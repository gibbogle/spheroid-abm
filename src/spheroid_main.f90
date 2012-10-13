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
integer :: jstep, hour, ntot, ncog, inflow
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

inbuflen = len(infile)
outbuflen = len(outfile)
write(*,*) 'call execute'
write(nfrun,*) 'infile: ',infile
write(nfrun,*) 'outfile: ',outfile
call execute(ncpu,infile,inbuflen,outfile,outbuflen)
!call cpu_time(t1)
write(*,*) 'did execute: nsteps: ',nsteps
do jstep = 1,Nsteps
!	write(*,*) 'jstep: ',jstep
	call simulate_step(res)
	if (res /= 0) then
		write(*,*) 'Error exit'
		stop
	endif
	if (mod(jstep,240) == 0) then
		call get_summary(summarydata)
!		hour = summaryData(1)
!		ntot = summaryData(3)
!		ncog = summaryData(5)
!		inflow = summaryData(7)
!		vasc = summaryData(8)/100.
!		write(*,'(4(a,i6),a,f6.2)') 'Hour: ',hour,' ncells: ',ntot,' ncog: ',ncog,' inflow/hr: ',inflow,' vasc: ',vasc	
!		write(nfrun,'(4(a,i6),a,f6.2)') 'Hour: ',hour,' ncells: ',ntot,' ncog: ',ncog,' inflow/hr: ',inflow,' vasc: ',vasc	
	endif
enddo
!call cpu_time(t2)
!write(*,*) 'time: ',t2-t1
end

