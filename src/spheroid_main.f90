
! Main program
!-----------------------------------------------------------------------------------------
PROGRAM spheroid_main
use spheroid_mod
use global
implicit none
integer :: ncpu, res, summarydata(100)
character*(128) :: infile, outfile, runfile
character*(64) :: travelfile = 'travel_time_dist.out'
integer :: status, nlen, cnt, i, inbuflen, outbuflen
integer :: jstep, hour, ntot, ncog, inflow, irun, i_hypoxia_cutoff,i_growth_cutoff, nsumm_interval
character*(128) :: b, c, progname
real :: vasc
real(8) :: t1, t2

runfile = 'spheroid_main.out'
open(nfrun,file=runfile,status='replace')
call disableTCP

outfile = 'spheroid_main.res'

call get_command (b, nlen, status)
if (status .ne. 0) then
    write (*,*) 'get_command failed with status = ', status
    stop
end if
!write (*,*) 'command line = ', b(1:len)
call get_command_argument (0, c, nlen, status)
if (status .ne. 0) then
    write (*,*) 'Getting command name failed with status = ', status
    stop
end if
!write (*,*) 'command name = ', c(1:len)
progname = c(1:nlen)
cnt = command_argument_count ()
!write (*,*) 'number of command arguments = ', cnt
if (cnt < 2) then
    write(*,*) 'Use: ',trim(progname),' num_cpu input_file'
    stop
endif

do i = 1, cnt
    call get_command_argument (i, c, nlen, status)
    if (status .ne. 0) then
        write (*,*) 'get_command_argument failed: status = ', status, ' arg = ', i
        stop
    end if
    if (i == 1) then
!        read(c(1:len),'(i)') ncpu
        read(c(1:nlen),*) ncpu															! --> ncpu
        write(*,*) 'Requested threads: ',ncpu
    elseif (i == 2) then
        infile = c(1:nlen)																! --> infile
        write(*,*) 'Input file: ',infile
    elseif (i == 3) then
        outfile = c(1:nlen)																! --> outfile
        write(*,*) 'Output file: ',outfile
    endif
end do

!call get_dimensions(NX,NY,NZ,nsteps,DELTA_T, MAX_CHEMO, cused);
i_hypoxia_cutoff = 3
i_growth_cutoff = 1
do irun = 1,1
	write(*,*) 'irun: ',irun
	inbuflen = len(infile)
	outbuflen = len(outfile)
	write(*,*) 'call execute'
	write(nfrun,*) 'infile: ',infile
	write(nfrun,*) 'outfile: ',outfile
	call execute(ncpu,infile,inbuflen,outfile,outbuflen)
	!call cpu_time(t1)
	t1 = wtime()
	write(*,*) 'did execute: nsteps, DELTA_T: ',nsteps, DELTA_T
	nsumm_interval = (60*60)/DELTA_T   ! number of time steps per hour
	write(*,*) 'nsumm_interval: ',nsumm_interval
	do jstep = 1,Nsteps
!		write(*,*) 'jstep: ',jstep
		call simulate_step(res)
		if (mod(jstep,nsumm_interval) == 0) then
			call get_summary(summarydata,i_hypoxia_cutoff,i_growth_cutoff)
		endif
		if (res /= 0) then
			write(*,*) 'Error exit: ',res
			stop
		endif
	enddo
	call terminate_run(res)
	!call cpu_time(t2)
	t2 = wtime()
	write(*,*) 'time: ',t2-t1
enddo
end

