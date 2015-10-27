#
# A ruby script to estimate memory requirement
#
# ruby mem.rb

def mem(nc, nnode, pm_factor= 3)
  # nc: number of particles per dim (np = nc^3)
  # nnode: number of MPI nodes
  # pm_factor: nc_particle_mesh = pm_factor * nc
  
  part_alloc_factor= 1.5

  puts "Memory estimate for nc= #{nc}, #{nnode} nodes" 
  sizep= 4*3*4 + 8 + 3*4   # particle and force
  sizep_min= 4*(3+3+2)   # position + velocity + id
  mbyte= 1024**2
  gbyte= 1024**3

  nx= nc/nnode
  if nc % nnode != 0 then
    puts "particle nc=#{nc} is not divisible by #{nnode}"
    nx= nc/nnode + 1
  end
  puts "particle nx= #{nx}"

  pm_nx= pm_factor*nc/nnode
  if pm_factor*nc % nnode != 0 then
    puts "PM nc=#{pm_factor*nc} is not divisible by #{nnode}"
    pm_nx= pm_factor*nc/nnode
  end
  puts "PM nx= #{pm_nx}"

  np= part_alloc_factor*nc*nc*(nx+1)         # number of particles per node
  mp= part_alloc_factor*nc*nc*(nx+1)*sizep   # RAM for particles
  mb= nc*nc*(1)*sizep;                       # particle communication buffer
  mm= (pm_factor*nc)**2*pm_nx*4              # PM mesh density -> force
  mm2= (pm_factor*nc)**2*(nx/2+1)*2*4        # PM mesh for density

  printf("%.1f Mbytes for %d particle\n", mp/mbyte, np)
  printf("%.1f Mbytes for communication buffer\n", mb/mbyte)
  printf("%.1f Mbytes for one PM mesh (mem1)\n", mm/mbyte)

  # memory for LPT
  m_lpt= nc**2*nx*4*12
  printf("%.1f Mbytes for 2LPT (mem1)\n", m_lpt/mbyte)
  
  # memory for snapshot
  m_snp= part_alloc_factor*nc*nc*(nx+1)*sizep_min

  printf("%.1f Mbytes for snapshot (mem2)\n", m_snp/mbyte)

  # memory shared for PM 2LPT/Particle Mesh/FoF halo finder/snapshot
  mem1= [m_lpt, mm].max
  mem2= [m_snp, mm2].max

  printf("%.1f Mbytes for mem1 (PM/LPT)\n", mem1/mbyte)
  printf("%.1f Mbytes for mem2 (PM/snapshot)\n", mem2/mbyte)

  printf("%.3f Gbytes per core in total (particle+buffer+mem1+mem2)\n", (mp + mb + mm + mm2)/gbyte)
  #printf("%.3f Gbyte for particle output\n", nc**3.0*32/gbyte) # long id
  printf("\n")
end

# WiZ-COLA configuration
# 1296^3 particles and 3 times finer grid

# 216 cores for g2 at Swinburne University (4GB RAM per core)
mem(1296, 216, 3)

# 432 cores necessary for 2GB RAM per core
mem(1296, 432, 3) 

