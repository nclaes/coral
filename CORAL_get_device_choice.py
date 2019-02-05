import numpy as np
import pyopencl as cl
from pyopencl import array
import sys
import subprocess
import time
import os
import argparse

def main():
	all_platforms = []
	platforms = cl.get_platforms()
	for pltm in platforms:
		if(pltm.name[:12] != 'Experimental'):	
			if(pltm.get_devices(cl.device_type.CPU) != []):
				all_platforms.append(pltm.get_devices(cl.device_type.CPU))

	for pltm in platforms:
		if(pltm.name[:12] != 'Experimental'):	
			if(pltm.get_devices(cl.device_type.GPU) != []):
				all_platforms.append(pltm.get_devices(cl.device_type.GPU))

	all_devices = []
	print('----------------------------------------------------------------------')
	print('Following is the list of platforms:')
	print('{:<15}'.format('PLATFORM NO.'),':','PLATFORM NAME/S\n')
	for i in range(len(all_platforms)):
		print('{:^15}'.format(i),':',all_platforms[i],'\n')
		for j in all_platforms[i]:
			# print(j)
			all_devices.append(j)
	print('----------------------------------------------------------------------\n')
	print('Following is the list of all devices:')
	print('{:<15}'.format('DEVICE NO.'),':','DEVICE NAME\n')
	for i in range(len(all_devices)):
		print('{:^15}'.format(i),':',all_devices[i],'\n')
	print('----------------------------------------------------------------------\n')
	# for i in range(len(all_platforms)):
	# 	if(all_devices[0] in all_platforms[i]):
	# 		print(all_devices[0],'in', all_platforms[i])


if __name__ == '__main__':
	main()
