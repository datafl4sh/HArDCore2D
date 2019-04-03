PROGRAM compute_rates

  IMPLICIT NONE

  INTEGER :: nbtests, i, io
  CHARACTER :: empty
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: meshsize, regmail, erl2p, erl2gradp, erlinfp, erflux, &
		ersharp, ersharpPe, Pemax, h1error
	CHARACTER (LEN=50) :: nomtest
	CHARACTER (LEN=400) :: commande
	
	! Rate file created by execseries.sh
	CHARACTER (LEN=*), PARAMETER :: rates_file="outputs/data_rates.dat"

	!
	! Read rate file
	!
	! nb de tests
	OPEN(10,file=rates_file,status="old")
	nbtests = -1
	DO
		READ(10,*,iostat=io) empty
		IF (io/=0 .OR. nbtests > 100) exit
		nbtests = nbtests+1
	ENDDO
	REWIND(10)

	! On lit les donnees
  ALLOCATE(meshsize(nbtests),regmail(nbtests),erl2p(nbtests),h1error(nbtests),erL2gradp(nbtests),erlinfp(nbtests),&
		erflux(nbtests),erSharp(nbtests),erSharpPe(nbtests),Pemax(nbtests))
	READ(10,*) empty
  DO i=1,nbtests
    READ(10,*) meshsize(i),erl2p(i),h1error(i)
		regmail(i)=1D0
  ENDDO
  CLOSE(10)

	!
	! Sortie a l'ecran
	!
	print *,"\n ---- Compute convergence rates ----"
	print "(A)","L2 error : h, error, rate, (mesh reg)"
	CALL output_rates(nbtests,meshsize,regmail,erl2p,20,"$L^2$ error on $p$")
	print "(A)","-----------------------------------------\nH1 error : h, error, rate, (mesh reg)"
	CALL output_rates(nbtests,meshsize,regmail,h1error,20,"H1 error on $p$")

	
	CLOSE(20)	

  DEALLOCATE(meshsize,erl2p,erl2gradp,erlinfp,erflux)

!---------------------------------------!
	CONTAINS

	SUBROUTINE output_rates(nbtests,meshsize,regmail,error,unitfich,type_error)
	!************************************!
	!
	!************************************!
	INTEGER :: nbtests, unitfich
	DOUBLE PRECISION, DIMENSION(:) :: meshsize, regmail, error
	CHARACTER (LEN=*) :: type_error

	DOUBLE PRECISION :: prev_meshsize, prev_er, rate


	prev_meshsize = 0D0
	prev_er = 0D0
	DO i=1,nbtests
		IF (prev_meshsize > 0D0) THEN
	    rate = (log(error(i))-log(prev_er))/(log(meshsize(i))-log(prev_meshsize))
			print "(ES9.2,2X,ES9.2,4X,F5.3,8X,A,ES8.2,A)",meshsize(i),error(i),rate,"(",regmail(i),")"

		ELSE
			print "(ES9.2,2X,ES9.2,17X,A,ES8.2,A)",meshsize(i),error(i),"(",regmail(i),")"

		ENDIF
		prev_meshsize = meshsize(i)
		prev_er = error(i)
	ENDDO

	END SUBROUTINE output_rates

END PROGRAM compute_rates
