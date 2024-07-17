This program, wemova, performs a WEighted MOVing Average on xyz
file. You need to give the name of the xyz file. The weights are
given in the kernel.dat file. The resulting file has a prefix
"filt_". Note, that the resultinf file is shirter than the
original.

In the example directory, you can find a linux binary, 
the kernel file and a short trajectory file. You can run 
the example with this command on linux:
"chmod +x wemova; printf water_traj.xyz | ./wemova"
