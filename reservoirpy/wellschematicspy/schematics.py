
# Source Code taken from:
#https://github.com/kinverarity1/well-schematics
import matplotlib.pyplot as plt
from matplotlib import patches as mpatches
from matplotlib import transforms as mtransforms
import numpy as np

PZONE_MAPPING = {
    "OH": "open hole",
    "open hole": "open hole",
    "open-hole": "open hole",
    "open": "open hole",
    "S": "wirewound screen",
    "screen": "wirewound screen",
    "wirewound screen": "wirewound screen",
    "SC": "slotted casing",
    "slots": "slotted casing",
    "slotted": "slotted casing",
    "slotted casing": "slotted casing",
}

class well_schema:
    def __init__(self, **kwargs):
        self.open_hole = kwargs.pop('open_hole',None)
        self.casing = kwargs.pop('casing',None)
        self.completion = kwargs.pop('completion',None)

    #Properties
    @property
    def open_hole(self):
        return self._open_hole

    @open_hole.setter
    def open_hole(self,value):
        if value is not None:
            assert isinstance(value,dict)
            for oh in value:
                assert isinstance(value[oh], dict)
                assert all(i in list(value[oh].keys()) for i in ['top','bottom','diameter'])
        self._open_hole = value


    @property
    def casing(self):
        return self._casing

    @casing.setter
    def casing(self,value):
        if value is not None:
            assert isinstance(value,dict)
            for oh in value:
                assert isinstance(value[oh], dict)
                assert all(i in list(value[oh].keys()) for i in ['type','top','bottom','diameter'])
                if 'cement' in list(value[oh].keys()):
                    assert isinstance(value[oh]['cement'],list)
                    for c in value[oh]['cement']:
                        assert isinstance(c,dict)
                        assert all(i in list(c.keys()) for i in ['top','bottom','oh'])

                if 'perforations' in list(value[oh].keys()):
                    assert isinstance(value[oh]['perforations'],list)
                    for c in value[oh]['perforations']:
                        assert isinstance(c,dict)
                        assert all(i in list(c.keys()) for i in ['top','bottom','oh'])
        self._casing = value

    @property
    def completion(self):
        return self._completion

    @completion.setter
    def completion(self,value):
        if value is not None:
            assert isinstance(value,dict)
            for oh in value:
                assert isinstance(value[oh], dict)
                assert all(i in list(value[oh].keys()) for i in ['type','top','bottom','diameter'])
        self._completion = value

    def max_diameter(self):
        d = []

        if self.open_hole is not None:
            d.extend([self.open_hole[o]['diameter'] for o in self.open_hole])
        if self.casing is not None:
            d.extend([self.casing[c]['diameter'] for c in self.casing])
        if self.completion is not None:
            d.extend([self.completion[c]['diameter'] for c in self.completion])
        
        return np.array(d).max()

    def unique_diameter(self):
        d = []

        if self.open_hole is not None:
            d.extend([self.open_hole[o]['diameter'] for o in self.open_hole])
        if self.casing is not None:
            d.extend([self.casing[c]['diameter'] for c in self.casing])
        if self.completion is not None:
            d.extend([self.completion[c]['diameter'] for c in self.completion])
        
        return np.unique(np.array(d))

    def top(self):
        d = []

        if self.open_hole is not None:
            d.extend([self.open_hole[o]['top'] for o in self.open_hole])
        if self.casing is not None:
            d.extend([self.casing[c]['top'] for c in self.casing])
        if self.completion is not None:
            d.extend([self.completion[c]['top'] for c in self.completion])
        
        return np.array(d).min()

    def bottom(self):
        d = []

        if self.open_hole is not None:
            d.extend([self.open_hole[o]['bottom'] for o in self.open_hole])
        if self.casing is not None:
            d.extend([self.casing[c]['bottom'] for c in self.casing])
        if self.completion is not None:
            d.extend([self.completion[c]['bottom'] for c in self.completion])
        
        return np.array(d).max()


    def plot(self,
        ax=None,
        tight_layout=True,
        dtick=True,
        xtick = True,
        lims=None,
        pipe_width=0.08,
        hatch_density=3,
        oh_kw={},
        csg_kw={}
    ):
        if ax is None:
            fig = plt.figure(figsize=(4, 9))
            ax = fig.add_subplot(111)

        t = mtransforms.blended_transform_factory(ax.transAxes, ax.transData)
        patches = []

        di_factor = self.max_diameter()

        def_oh_kw = {
            'color': '#cfd4d3',
            'fill':True,
            'hatch':None
        }    

        for (k,v) in def_oh_kw.items():
            if k not in oh_kw:
                oh_kw[k]=v

        def_csg_kw = {
            'color': 'k',
            'pipe_width':0.08,
            'shoe_scale':5
        }    

        for (k,v) in def_csg_kw.items():
            if k not in csg_kw:
                csg_kw[k]=v

        #Open Hole
        if self.open_hole is not None:
            for o in self.open_hole:
                top = self.open_hole[o]['top']
                bottom = self.open_hole[o]['bottom']
                length = bottom - top
                diameter = self.open_hole[o]['diameter']
                color = self.open_hole[o].pop('color','#cfd4d3')
                hatch = self.open_hole[o].pop('hatch',None)
                oh_patch = mpatches.Rectangle(
                    (0.5*(1-diameter/di_factor), top),
                    (0.5*(1+diameter/di_factor)) - (0.5*(1-diameter/di_factor)),
                    length,
                    facecolor=color,
                    transform=t,
                    hatch = hatch
                )
                patches.append(oh_patch)

        if self.casing is not None:
            for c in self.casing:
                ctype = self.casing[c]['type']
                top = self.casing[c]['top']
                bottom = self.casing[c]['bottom']
                length = bottom - top
                diameter = self.casing[c]['diameter']
                pipe_width=self.casing[c].pop('pipe_width',0.04)
                shoe_scale=self.casing[c].pop('shoe_scale',5)
                color = self.casing[c].pop('color','k')

                xl =  0.5*(1-diameter/di_factor)
                xr =  0.5*(1+diameter/di_factor)
                seg_left = mpatches.Rectangle(
                    (xl, top), 
                    pipe_width, 
                    length, 
                    facecolor=color,
                    transform=t,
                )
                seg_right = mpatches.Rectangle(
                    (xr-pipe_width, top),
                    pipe_width,
                    length,
                    facecolor=color,
                    transform=t,
                )
                
                #Shoe
                left_shoe = np.array([[xl,0],[xl,-1*shoe_scale],[xl-pipe_width,0]])
                left_shoe[:,1] = left_shoe[:,1] + bottom
                right_shoe = np.array([[xr,0],[xr,-1*shoe_scale],[xr+pipe_width,0]])
                right_shoe[:,1] = right_shoe[:,1] + bottom
                ls = mpatches.Polygon(left_shoe, color='k')
                rs = mpatches.Polygon(right_shoe, color='k')

                patches.extend([seg_left,seg_right,ls,rs])

                if 'cement' in self.casing[c]:
                    for cem in self.casing[c]['cement']:
                        cement_color = cem.pop('color','#adadad')
                        cement_hatch = cem.pop('hatch','/')
                        cement_oh = cem['oh']
                        cement_top = cem['top']
                        cement_bottom = cem['bottom']

                        cl =  0.5*(1-cement_oh/di_factor)
                        cement_width = xl - cl
                        cement_length = cement_bottom - cement_top

                        cem_left = mpatches.Rectangle(
                            (cl, cement_top), 
                            cement_width, 
                            cement_length, 
                            facecolor=cement_color,
                            transform=t,
                            hatch = cement_hatch
                        )

                        cem_right = mpatches.Rectangle(
                            (xr, cement_top), 
                            cement_width, 
                            cement_length, 
                            facecolor=cement_color,
                            transform=t,
                            hatch = cement_hatch
                        )
                        patches.extend([cem_left,cem_right])
                
                if 'perforations' in self.casing[c]:
                    for perf in self.casing[c]['perforations']:
                        perf_color = perf.pop('color','#030302')
                        perf_hatch = perf.pop('hatch','*')
                        perf_scale = perf.pop('scale',1)
                        perf_penetrate = perf.pop('penetrate',1.05)
                        perf_oh = perf['oh']
                        perf_top = perf['top']
                        perf_bottom = perf['bottom']

                        pl =  0.5*(1-perf_oh*perf_penetrate/di_factor)
                        pr =  0.5*(1+perf_oh*perf_penetrate/di_factor)
                        
                        for i in np.arange(perf_top,perf_bottom,perf_scale):
                            left_perf = np.array([[pl,perf_scale/2],[xl,perf_scale],[xl,0]])
                            left_perf[:,1] = left_perf[:,1] + i
                            right_perf = np.array([[pr,perf_scale/2],[xr,perf_scale],[xr,0]])
                            right_perf[:,1] = right_perf[:,1] + i

                            lp = mpatches.Polygon(
                                left_perf, 
                                color=perf_color,
                                hatch=perf_hatch
                            )
                            rp = mpatches.Polygon(
                                right_perf, 
                                color=perf_color,
                                hatch=perf_hatch
                            )
                            patches.extend([lp,rp])

        for patch in patches:
            ax.add_artist(patch)

        ax.grid(False)
        for side in ["left", "right", "bottom", "top"]:
            ax.spines[side].set_visible(True)
        if not dtick:
            ax.yaxis.set_ticks_position("none")
        #ax.set_facecolor("white")
        if xtick:
            di = self.unique_diameter()
            ax.set_xticks(np.concatenate([0.5*(1-di/di_factor),0.5*(1+di/di_factor)]))
            ax.set_xticklabels(np.around(np.concatenate([di,di]),decimals=1))
        else:
            ax.set_xticklabels([])
        ax.xaxis.set_label_position("top")
        ax.xaxis.tick_top()
        ax.set_xlim([0,1])

        if lims is None:
            ax.set_ylim([self.bottom(),self.top()])
        else:
            ax.set_ylim([lims[1],lims[0]])

        if tight_layout:
            ax.figure.tight_layout()

        return patches





def plot_schematic(
    segments,
    ax=None,
    tight_layout=True,
    depth_tick_markers=False,
    pipe_width=0.08,
    hatch_density=3,
):
    """Draw casing in a well which is a single diameter construction.
    Args:
        segments (sequence of dicts): each dict should be in the
            form ``{"type": <str>, "top": <float>, "bottom": <float>}``.
            The "type" should be either "casing", "pipe", "blank", or "sump",
            or a production zone type (either "screen", "slotted casing" or
            "open hole"). "top" and "bottom" are the top and bottom of each
            segment.
        ax (matplotlib.Axes): to draw in
        tight_layout (bool): run tight_layout() on ax.figure to rearrange
            things to fit.
        depth_tick_markers (bool): show tick markers for the vertical
            depth axis. Labels will always appear.
        pipe_width (float): width of pipe
        hatch_density (int): density of screen hatching
    Returns: a list of the artists created.
    """
    if ax is None:
        fig = plt.figure(figsize=(4, 9))
        ax = fig.add_subplot(111)

    t = mtransforms.blended_transform_factory(ax.transAxes, ax.transData)
    patches = []

    di = np.array([s['diameter'] for s in segments])
    di_factor = di.max()
    print(di_factor)
    for segment in segments:
        seg_type = segment["type"]
        seg_from = segment["top"]
        seg_to = segment["bottom"]
        seg_dia = segment['diameter']
        seg_length = seg_to - seg_from

        if seg_type in PZONE_MAPPING:
            seg_type = PZONE_MAPPING[seg_type]
            if seg_type == "wirewound screen":
                hatch = "-" * hatch_density
            elif seg_type == "slotted casing":
                hatch = "/" * hatch_density
            else:
                hatch = None

            seg_left = mpatches.Rectangle(
                (0.5*(1-seg_dia/di_factor), seg_from),
                pipe_width * 0.9,
                seg_length,
                facecolor="k",
                fill=False,
                hatch=hatch,
                transform=t,
            )
            seg_right = mpatches.Rectangle(
                (0.5*(1+seg_dia/di_factor)-pipe_width, seg_from),
                pipe_width * 0.9,
                seg_length,
                facecolor="k",
                fill=False,
                hatch=hatch,
                transform=t,
            )
        else:
            seg_type = "pipe"
            seg_left = mpatches.Rectangle(
                (0.5*(1-seg_dia/di_factor), seg_from), 
                pipe_width, 
                seg_length, 
                facecolor="k", 
                transform=t
            )
            seg_right = mpatches.Rectangle(
                (0.5*(1+seg_dia/di_factor)-pipe_width, seg_from),
                pipe_width,
                seg_length,
                facecolor="k",
                transform=t,
            )
        patches += [seg_left, seg_right]

    for patch in patches:
        ax.add_artist(patch)

    ax.grid(False)
    for side in ["left", "right", "bottom", "top"]:
        ax.spines[side].set_visible(True)
    if not depth_tick_markers:
        ax.yaxis.set_ticks_position("none")
    ax.set_facecolor("white")
    
    ax.set_xticks(np.concatenate([0.5*(1-di/di_factor),0.5*(1+di/di_factor)]))
    ax.set_xticklabels(np.around(np.concatenate([di,di]),decimals=1))
    ax.xaxis.set_label_position("top")
    ax.xaxis.tick_top()
    ax.set_xlim([-0.2,1.2])
    ax.set_ylim(
        max([s["bottom"] for s in segments]) + 1, min([s["top"] for s in segments]) - 1
    )
    if tight_layout:
        ax.figure.tight_layout()

    return patches