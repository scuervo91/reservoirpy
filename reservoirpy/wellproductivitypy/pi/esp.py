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
	head_coef: List[float] = Field(..., min_items=6,max_items=6)
	power_coef: Optional[List[float]] = Field(...,min_items=6,max_items=6)
	ref_stage : int = Field(1, ge=1)
	ref_freq : float = Field(60, ge=0)
	n : int = Field(20, ge = 10)


	def performance(self, stages:list[int]=None, freq:List[float]=None, viscosity_correction:bool=False):

		if stages is None:
			stages = np.array(self.ref_stage)
		else:
			stages = np.atleast_1d(stages)

		if freq is None:
			freq = np.atleast_2d(self.ref_freq).T
		else:
			freq = np.atleast_2d(freq).T 


		# Capacity Ranges ref and adjusted by frecuency
		capacity_range_ref = np.linspance(0,self.aof,self.n)
		capacity_range = capacity_range_ref * (freq/self.ref_freq)

		#Polynomials Head
		head_poly = np.poly1d(self.head_coef)

		head_ref = head_poly(capacity_range_ref)
		head = head_ref * stages * np.power(freq/self.ref_freq,2)

		# Polynomials Power

		power_poly = np.poly1d(self.power_poly)

		power_ref = power_poly(capacity_range_ref)
		power = power_ref * stages * np.power(freq/self.ref_freq,3)

		hydro_power = capacity_range_ref * head_ref * 62.4 * 1.1816e-7
		efficiency = hydro_power / power_ref


		#bounds
		cap_min = self.min_flow * (freq/self.ref_freq)
		cap_max = self.max_flow * (freq/self.ref_freq)

		head_min = head_poly(self.min_flow) * stages * np.power(freq/self.ref_freq,2)
		head_min = head_poly(self.max_flow) * stages * np.power(freq/self.ref_freq,2)

		return 














