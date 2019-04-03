PROGRAM calcul_taux

  IMPLICIT NONE

  INTEGER :: nbtests, i, io
  CHARACTER :: vide
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: pas, regmail, erl2p, erl2gradp, erlinfp, erflux, &
		ersharp, ersharpPe, Pemax, h1error
	CHARACTER (LEN=50) :: nomtest
	CHARACTER (LEN=400) :: commande
	
	! Fichier de taux cree par execseries.sh
	CHARACTER (LEN=*), PARAMETER :: fichier_taux="outputs/data_taux.dat"

	!
	! On commence par lire le fichier de taux
	!
	! nb de tests
	OPEN(10,file=fichier_taux,status="old")
	nbtests = -1
	DO
		READ(10,*,iostat=io) vide
		IF (io/=0 .OR. nbtests > 100) exit
		nbtests = nbtests+1
	ENDDO
	REWIND(10)

	! On lit les donnees
  ALLOCATE(pas(nbtests),regmail(nbtests),erl2p(nbtests),h1error(nbtests),erL2gradp(nbtests),erlinfp(nbtests),&
		erflux(nbtests),erSharp(nbtests),erSharpPe(nbtests),Pemax(nbtests))
	READ(10,*) vide
  DO i=1,nbtests
!    READ(10,*) pas(i),regmail(i),erl2p(i),erl2gradp(i),erlinfp(i),erflux(i),erSharp(i),&
!			erSharpPe(i),Pemax(i),nomtest
    READ(10,*) pas(i),erl2p(i),h1error(i)
		regmail(i)=1D0
  ENDDO
  CLOSE(10)

	!
	! Sortie a l'ecran
	!
	print *,"\n ---- Calcul des taux de convergence ----"
	! Taux de cv L2 sur p
	print "(A)","L2 error : h, error, rate, (mesh reg)"
	CALL sortie_taux(nbtests,pas,regmail,erl2p,20,"$L^2$ error on $p$")
	print "(A)","-----------------------------------------\nH1 error : h, error, rate, (mesh reg)"
	CALL sortie_taux(nbtests,pas,regmail,h1error,20,"H1 error on $p$")
	! Taux de cv L2 sur gradp
!	print "(A)","Erreur L2 gradp : pas, erreur, taux, (reg maillage)"
!	CALL sortie_taux(nbtests,pas,regmail,erl2gradp,20,"$L^2$ error on $\\nabla p$")
!	! Taux de cv Linf sur p
!	print "(A)","Erreur Linf p : pas, erreur, taux, (reg maillage)"
!	CALL sortie_taux(nbtests,pas,regmail,erlinfp,20,"$L^\\infty$ error on $p$")
!	! Taux de cv sur le flux
!	print "(A)","Erreur sur le flux : pas, erreur, taux, (reg maillage)"
!	CALL sortie_taux(nbtests,pas,regmail,erflux,20,"Flux error")
!	! Taux de cv sur la norme sharp
!	print "(A)","Erreur sur norme sharp : pas, erreur, taux, (reg maillage)"
!	CALL sortie_taux(nbtests,pas,regmail,erSharp,20,"Flux error")
!	! Taux de cv sur la norme sharp / min(1,Pe^{1/2})
!	print "(A)","Erreur sur norme sharp / min(1,Pe^{1/2}) : pas, erreur, taux, (reg maillage)"
!	CALL sortie_taux(nbtests,pas,regmail,erSharpPe,20,"Flux error")

	
	CLOSE(20)	

  DEALLOCATE(pas,erl2p,erl2gradp,erlinfp,erflux)

!---------------------------------------!
	CONTAINS

	SUBROUTINE sortie_taux(nbtests,pas,regmail,erreur,unitfich,type_erreur)
	!************************************!
	!	Affiche a l'ecran et dans le fichier latex
	!	la liste des pas, erreur et taux correspondants
	!
	!************************************!
	INTEGER :: nbtests, unitfich
	DOUBLE PRECISION, DIMENSION(:) :: pas, regmail, erreur
	CHARACTER (LEN=*) :: type_erreur

	DOUBLE PRECISION :: pas_prec, er_prec, taux


	pas_prec = 0D0
	er_prec = 0D0
	DO i=1,nbtests
		IF (pas_prec > 0D0) THEN
	    taux = (log(erreur(i))-log(er_prec))/(log(pas(i))-log(pas_prec))
			print "(ES9.2,2X,ES9.2,4X,F5.3,8X,A,ES8.2,A)",pas(i),erreur(i),taux,"(",regmail(i),")"

		ELSE
			print "(ES9.2,2X,ES9.2,17X,A,ES8.2,A)",pas(i),erreur(i),"(",regmail(i),")"

		ENDIF
		pas_prec = pas(i)
		er_prec = erreur(i)
	ENDDO

	END SUBROUTINE sortie_taux

END PROGRAM calcul_taux
