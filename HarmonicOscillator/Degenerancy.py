#!/usr/local/opt/python/bin/python3.7
__author__ = '@Tssp'

# E_n = hw(N + 3/2) with N = nx + ny + nz. Degenerancy for every combination with same N
def get_states(N):
	out = list()
	for nx in range(N+1):
		for ny in range(N+1):
			for nz in range(N+1):
				if (nx + ny + nz) == N:
					out.append("(" + str(nx) + ", " + str(ny) + ", " + str(nz) + ")")
				else:
					continue
	fraction = (N+3/2).as_integer_ratio()
	print(f"E = {fraction[0]}/{fraction[1]}hw {len(out)} states:\n {out}\n")
	return out

for N in range(4):
	get_states(N)
