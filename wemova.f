c      This program, wemova, performs a WEighted MOVing Average on xyz
c      file. You need to give the name of the xyz file. The weights are
c      given in the kernel.dat file. The resulting file has a prefix
c      "filt_". Note, that the resultinf file is shirter than the
c      original.
c      email: madarasz.adam@ttk.hu

      program wemova

      implicit none

      integer ndec,ablak,nblock,natoms
      integer i,j,k,m,ipos,tav,midpos
      integer ierror,iline
      DOUBLE PRECISION kernel(0:100000),PI,temp
      DOUBLE PRECISION sum,timestep
      DOUBLE PRECISION, allocatable :: xtot(:,:),ytot(:,:),ztot(:,:)
      DOUBLE PRECISION, allocatable :: cx(:),cy(:),cz(:)
      character*120 xyzfile,outxyzfile,stratoms
      character(300), allocatable :: title(:)  
      character(5), allocatable :: nam(:)  
      LOGICAL :: kernel_file_exist
      character*120 kernel_file_name

      kernel_file_name="kernel.dat"

      PI=4.D0*DATAN(1.D0)

c
c     get the base name of user specified input structures
c
      write (*,10)
  10    format (/,' Name of the input xyz file: ',$)
      read (*,20)  xyzfile
   20    format (a120)

      write (*,*) xyzfile

c
c     get the number of frames
c

      OPEN (50, file = xyzfile)
      read(50,'(i10)',iostat=ierror) natoms

      if ( ierror .ne. 0 ) then

            write (*,*) "ERROR: Incorrect xyz file format. ",
     +          "The first line is not an integer ",
     +          "that should mean the number of atoms. Program stops."
            stop

      endif

      i = 1
      DO
       READ (50,*, END=11)
       i = i + 1
      END DO
   11 CLOSE (50)

      if (mod(i,(natoms+2)) .ne. 0) then

            write (*,*) "WARNING: The number of lines does not ",
     +          "match the number of atoms."
            stop

      endif 

      nblock=i/(natoms+2)

      write(*,*) 'Number of atoms: ',natoms
      write(*,*) 'Number of input frames: ',nblock

c      write (*,80)
c   80 format (/,' Timestep between frames',
c     &              ' in picoseconds:  ',$)
c   90 format (f20.0)
c      read (*,90) timestep
c      write (*,*) timestep

c Reading the kernel function

      INQUIRE(FILE=kernel_file_name, EXIST=kernel_file_exist)

      if (kernel_file_exist) then

        OPEN (51, file = kernel_file_name)

         j = 0
         DO
           READ (51,*, END=14) kernel(j)
           j = j + 1
         END DO
   14 CLOSE (51)

      else

          write(*,*) "kernel.dat file is not found. Program stops."
          stop

      endif

      ndec=j-1

      ablak=2*ndec+1

c Check the integral of the kernel function

      sum=kernel(0)/2.0d0

      do j=1,ndec

          sum=sum+kernel(j)

      end do

      sum=sum*2.0d0

      write(*,*) "The integral of the kernel function is: ",sum
      write(*,*) "before it is normalized to 1."
c Normalize the kernel function to 1.0

      do j=0,ndec

          write(*,*) kernel(j)

          kernel(j)=kernel(j)/sum

      end do

      write(*,*) 'Starting filtration'
c
c     perform dynamic allocation of some local arrays
c
      allocate (title(nblock))
      allocate (nam(natoms))
      allocate (cx(natoms))
      allocate (cy(natoms))
      allocate (cz(natoms))

      allocate (xtot(natoms,ablak))
      allocate (ytot(natoms,ablak))
      allocate (ztot(natoms,ablak))

      open (unit=50,file=xyzfile,status='old',action='read')  

      outxyzfile='filt_'//xyzfile

      open (unit=60,file=outxyzfile,status='unknown',action='write')  

      write(*,*) 'Name of the outputfile: ', outxyzfile
      write(*,*) 'Number of filtered frames: ', nblock-ablak+1

      iline=0

c
c     cycle over all pairs of snapshot frame blocks
c
      do i = 1, nblock

         read(50,'(a)') stratoms
         read(50,15) title(i)
   15    format (a300)

         iline=iline+2

         do j=1,natoms
            iline=iline+1
            read(50,*,iostat=ierror) nam(j),cx(j),cy(j),cz(j)  

            if ( ierror .ne. 0 ) then

                write (*,*) "ERROR: Problem reading file ",xyzfile
                write (*,*) "in line ",iline
                write (*,*)  "Program stops."
            stop

            endif

         end do  
 
         ipos=mod(i-1,ablak)+1
         midpos=mod(i-1-ndec,ablak)+1
         do j = 1, natoms
            xtot(j,ipos)=cx(j)
            ytot(j,ipos)=cy(j)
            ztot(j,ipos)=cz(j)
         end do

         if (i .ge. ablak) then

            do m=1, natoms
               cx(m)=0.0d0
               cy(m)=0.0d0
               cz(m)=0.0d0
            end do
               
            do k=1,ablak
               tav=min(abs(midpos-k),ablak+k-midpos,ablak+midpos-k)
               do m=1, natoms
                  cx(m)=cx(m)+xtot(m,k)*kernel(tav)
                  cy(m)=cy(m)+ytot(m,k)*kernel(tav)
                  cz(m)=cz(m)+ztot(m,k)*kernel(tav)
               end do
            end do

c
c     write output
c

      write(*,*) 'Writing frame', i-2*ndec
      
            write(60,'(a)') trim(stratoms)

            write(60,'(a)') trim(title(i-ndec))

            do k=1,natoms
               write(60,25) nam(k),cx(k),cy(k),cz(k)  
   25    format (a3,3f11.5)
            end do  

         end if
         
      end do

c
c     perform deallocation of some local arrays
c
      deallocate (title)
      deallocate (nam)
      deallocate (cx)
      deallocate (cy)
      deallocate (cz)

      deallocate (xtot)
      deallocate (ytot)
      deallocate (ztot)


      close(50)
      close(60)

      write(*,*) "wemova terminated normally."

      end program wemova

