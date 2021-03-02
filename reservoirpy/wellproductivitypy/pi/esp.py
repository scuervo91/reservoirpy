import numpy as np 
import pandas as pd 
from pydantic import BaseModel, Field
from typing import Optional, List


class Esp(BaseModel):
	name : str = Field(...)
	series : Optional[int] = Field(None, gt=0)
	model : Optional[str] = Field(None)
	min_csg : Optional[float] = Field(None)
	min_flow: float = Field(..., gt=0)
	max_flow: float = Field(..., gt=0)
	bep: float = Field(..., gt=0)
	aof : float = Field(..., gt=0)
	head_coef: List[float] = Field(..., min_items=3,max_items=7)
	power_coef: Optional[List[float]] = Field(...,min_items=3,max_items=7)
	ref_stage : int = Field(1, ge=1)
	ref_freq : float = Field(60, ge=0)
	n : int = Field(20, gt = 1)


	def head(self, stages:List[int]=None, freq:List[float]=None, viscosity_correction:bool=False):

		if stages is None:
			stages = np.array(self.ref_stage)
		else:
			stages = np.atleast_1d(stages)

		if freq is None:
			freq = np.atleast_2d(self.ref_freq).T
		else:
			freq = np.atleast_2d(freq).T 


		# Capacity Ranges ref and adjusted by frecuency
		capacity_range_ref = np.linspace(0,self.aof,self.n)
		capacity_range = capacity_range_ref * (freq/self.ref_freq)

		#Polynomials Head
		head_poly = np.poly1d(self.head_coef)

		head_ref = head_poly(capacity_range_ref)
		head = head_ref * stages * np.power(freq/self.ref_freq,2)

		freq_array = np.repeat(freq.reshape(-1),self.n)
		stage_array = np.full(freq_array.shape, stages)
		head_array = head.flatten()


		df = pd.DataFrame(
			np.column_stack((head_array,freq_array,stage_array)), 
			index = pd.Index(capacity_range.flatten(),name='capacity'),
			columns = ['head','frecuency','stages']
		)

		return df

	def operation_range(self, stages:List[int]=None, freq:List[float]=None, viscosity_correction:bool=False):

		if stages is None:
			stages = np.array(self.ref_stage)
		else:
			stages = np.atleast_1d(stages)

		if freq is None:
			freq = np.atleast_2d(self.ref_freq).T
		else:
			freq = np.atleast_2d(freq).T 

		head_poly = np.poly1d(self.head_coef)

		#bounds
		cap_min = self.min_flow * (freq/self.ref_freq)
		cap_max = self.max_flow * (freq/self.ref_freq)

		head_min = head_poly(self.min_flow) * stages * np.power(freq/self.ref_freq,2)
		head_max = head_poly(self.max_flow) * stages * np.power(freq/self.ref_freq,2)

		df = pd.DataFrame({
			'capacity': np.concatenate((cap_min,cap_max),axis=None),
			'head': np.concatenate((head_min,head_max), axis=None),
			'bound': np.repeat(['l','u'],freq.shape[0])
		})

		return df

	def power(self, stages:List[int]=None, freq:List[float]=None, viscosity_correction:bool=False):

		if stages is None:
			stages = np.array(self.ref_stage)
		else:
			stages = np.atleast_1d(stages)

		if freq is None:
			freq = np.atleast_2d(self.ref_freq).T
		else:
			freq = np.atleast_2d(freq).T 

		# Capacity Ranges ref and adjusted by frecuency
		capacity_range_ref = np.linspace(0,self.aof,self.n)
		capacity_range = capacity_range_ref * (freq/self.ref_freq)

		# Polynomials Power

		head_poly = np.poly1d(self.head_coef)

		head_ref = head_poly(capacity_range_ref)

		power_poly = np.poly1d(self.power_coef)

		power_ref = power_poly(capacity_range_ref)
		power = power_ref * stages * np.power(freq/self.ref_freq,3)

		hydro_power = capacity_range_ref * head_ref * 62.4 * 1.1816e-7
		efficiency = hydro_power / power_ref

		freq_array = np.repeat(freq.reshape(-1),self.n)
		power_array = power.flatten()
		efficiency_array = np.tile(efficiency.flatten(),freq.shape[0])
		stage_array = np.full(freq_array.shape, stages)

		df = pd.DataFrame(
			np.column_stack((power_array,efficiency_array,freq_array,stage_array)), 
			index = pd.Index(capacity_range.flatten(),name='capacity'),
			columns = ['power','efficiency','frecuency','stages']
		)

		return df



