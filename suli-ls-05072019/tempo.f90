	!Subrotina para calcular o 


	!!! Implementação 03/07/2017
	! Leonardo Romero Monteiro

	!!! Modificações
	! Leonardo Romero Monteiro - 103/07/2017

	SUBROUTINE tempo()

	USE vartempo
	USE velpre

	IMPLICIT NONE

	!===================================================================================================================

	!contadores
	integer :: i, j, k
	real :: aux1

	if ((t_tempo == 0) .or. ((t_tempo == 3) .and. (it == 1))) then ! Euler e primeiro do AB2
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx1
		u(i,j,k) = u(i,j,k) + dt*Fu(i,j,k)

	enddo
	enddo	
	enddo

	do k = 1, nz
	do j = 1, ny1
	do i = 1, nx
		v(i,j,k) = v(i,j,k) + dt*Fv(i,j,k)

	enddo
	enddo
	enddo

	do k = 1, nz1
	do j = 1, ny
	do i = 1, nx
		w(i,j,k) = w(i,j,k) + dt*Fw(i,j,k)

	enddo
	enddo
	enddo
		if (t_tempo == 3) then ! só para AB2
		fu0 = Fu
		fv0 = Fv
		fw0 = Fw
		endif

	elseif (t_tempo == 1) then  !RK 2

		if (tt == 1) then
		aux1 = dt ! já está modificado para ser 0.5 do dt original

			do k = 1, nz
			do j = 1, ny
			do i = 1, nx1
				u0(i,j,k)  = u(i,j,k)
				u(i,j,k)   = u(i,j,k) + aux1*Fu(i,j,k)
			enddo
			enddo	
			enddo

			do k = 1, nz
			do j = 1, ny1
			do i = 1, nx
				v0(i,j,k)  = v(i,j,k)
				v(i,j,k)   = v(i,j,k) + aux1*Fv(i,j,k)
			enddo
			enddo
			enddo

			do k = 1, nz1
			do j = 1, ny
			do i = 1, nx
				w0(i,j,k)  = w(i,j,k)
				w(i,j,k)   = w(i,j,k) + aux1*Fw(i,j,k)
			enddo
			enddo
			enddo

		elseif (tt == 2) then
		aux1 = dt

			do k = 1, nz
			do j = 1, ny
			do i = 1, nx1
				u(i,j,k) = u0(i,j,k) +aux1 * Fu(i,j,k)	
			enddo
			enddo	
			enddo

			do k = 1, nz
			do j = 1, ny1
			do i = 1, nx
				v(i,j,k) = v0(i,j,k) +aux1 * Fv(i,j,k)
			enddo
			enddo
			enddo

			do k = 1, nz1
			do j = 1, ny
			do i = 1, nx
				w(i,j,k) = w0(i,j,k) +aux1 * Fw(i,j,k)
			enddo
			enddo
			enddo

		endif


	elseif (t_tempo == 2) then  !RK 3

		if (tt == 1) then
		aux1 = dt ! já está modificado para ser 0.5 do dt original

			do k = 1, nz
			do j = 1, ny
			do i = 1, nx1
				u0(i,j,k)  = u(i,j,k)
				fu0(i,j,k) = Fu(i,j,k)		
				u(i,j,k)   = u(i,j,k) + aux1*Fu(i,j,k)
			enddo
			enddo	
			enddo

			do k = 1, nz
			do j = 1, ny1
			do i = 1, nx
				v0(i,j,k)  = v(i,j,k)
				fv0(i,j,k) = Fv(i,j,k)	
				v(i,j,k)   = v(i,j,k) + aux1*Fv(i,j,k)
			enddo
			enddo
			enddo

			do k = 1, nz1
			do j = 1, ny
			do i = 1, nx
				w0(i,j,k)  = w(i,j,k)
				fw0(i,j,k) = Fw(i,j,k)	
				w(i,j,k)   = w(i,j,k) + aux1*Fw(i,j,k)
			enddo
			enddo
			enddo

		elseif (tt == 2) then
		aux1 = 2.*dt

			do k = 1, nz
			do j = 1, ny
			do i = 1, nx1
				fu1(i,j,k) = Fu(i,j,k)		
				u(i,j,k)   = u0(i,j,k) - dt*fu0(i,j,k) + aux1*Fu(i,j,k)	
			enddo
			enddo	
			enddo

			do k = 1, nz
			do j = 1, ny1
			do i = 1, nx
				fv1(i,j,k) = Fv(i,j,k)	
				v(i,j,k)   = v0(i,j,k) - dt*fv0(i,j,k) + aux1*Fv(i,j,k)
			enddo
			enddo
			enddo

			do k = 1, nz1
			do j = 1, ny
			do i = 1, nx
				fw1(i,j,k) = Fw(i,j,k)	
				w(i,j,k)   = w0(i,j,k) - dt*fw0(i,j,k) + aux1*Fw(i,j,k)
			enddo
			enddo
			enddo

		elseif (tt == 3) then
		aux1 = dt/6.

			do k = 1, nz
			do j = 1, ny
			do i = 1, nx1
				u(i,j,k) = u0(i,j,k) +aux1 * (fu0(i,j,k) + 4.*fu1(i,j,k) + Fu(i,j,k)) 	
			enddo
			enddo	
			enddo

			do k = 1, nz
			do j = 1, ny1
			do i = 1, nx
				v(i,j,k) = v0(i,j,k) +aux1 * (fv0(i,j,k) + 4.*fv1(i,j,k) + Fv(i,j,k)) 	
			enddo
			enddo
			enddo

			do k = 1, nz1
			do j = 1, ny
			do i = 1, nx
				w(i,j,k) = w0(i,j,k) +aux1 * (fw0(i,j,k) + 4.*fw1(i,j,k) + Fw(i,j,k)) 	
			enddo
			enddo
			enddo


		endif




	elseif ( (t_tempo == 3) .and. (it .ne. 1)) then ! Euler e primeiro do AB2
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx1
		u(i,j,k) = u(i,j,k) + dt*0.5*(3*Fu(i,j,k)-fu0(i,j,k))
		fu0(i,j,k) = Fu(i,j,k)	
	enddo
	enddo	
	enddo

	do k = 1, nz
	do j = 1, ny1
	do i = 1, nx
		v(i,j,k) = v(i,j,k) + dt*0.5*(3*Fv(i,j,k)-fv0(i,j,k))
		fv0(i,j,k) = Fv(i,j,k)	
	enddo
	enddo
	enddo

	do k = 1, nz1
	do j = 1, ny
	do i = 1, nx
		w(i,j,k) = w(i,j,k) + dt*0.5*(3*Fw(i,j,k)-fw0(i,j,k))
		fw0(i,j,k) = Fw(i,j,k)	
	enddo
	enddo
	enddo



	endif


!==================================================================================================================
END SUBROUTINE tempo





