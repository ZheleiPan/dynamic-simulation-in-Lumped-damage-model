c
c*****************************************************************
c
c       Elemento finito con variables internas para MDC
c            Modelo Con Corrosi¢n
c
c       MDC
c       VERSION 1.0
c
c******************************************************************
c
c       SUBRUTINA DE INTERFAZ CON ABAQUS
c       VERSION 1.0
c       ESCRITA POR Scarlet, Ricardo, Carlos y Julio Flrez Lpez
c       Mayo 2021
c
c******************************************************************** 
c              
        subroutine uel(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,
     1              props,nprops,coords,mcrd,nnode,u,du,v,a,jtype,
     2              time,dtime,kstep,kinc,jelem,params,
     3              ndload,jdltyp,adlmag,predef,npredf,lflags,
     4              mlvarx,ddlmag,mdload,pnewdt,jprops,njprop,period)
c
c
c       IMPLICIT NONE
c*************************************************************
c              DECLARACION DE LAS VARIABLES ABAQUS
c************************************************************
c
c       para el significado de estas variables ver el manual
c       ABAQUS  versin 6.1 seccin 23.2.19
c
c       Las caracterçsticas del elemento y de la integracin numârica
c       son transmitidas por abaqus a pdp_local a travâs de la matriz
c
c       PAR¬METROS DE IDENTIFICACI‘N DEL ELEMENTO
c
c       jprops(3) nc =Nîmero de coordenadas que definen el nodo = 2
c       jprops(4) nfi=Nîmero de grados de libertad del nodo = 3
c       jprops(5) nv =Nîmero de variables internas con funcin de fluencia = 6
c       jprops(6) na =Nîmero de variables internas sin funcin de fluencia = 4
c       jprops(7) ny =Nîmero de fuerzas termodinamicas = 6
c
c       PAR¬METROS DE INTEGRACI‘N NUM–RICA
c       jprops(1) max_ite=maximo No.de iteraciones = 5
c       jprops(2) max_paso=maximo No.de pasos de integracion local = 1
c       props(1:nfi+nv+na)  Ceros numâricos de {R},{T},{G} = 13
c       params(1:3) Coeficientes del mâtodo de newmark
c
c       PROPIEDADES DEL ELEMENTO
c       props(1:nfi+nv+na)  propiedades del elemento = 13
c
c       VARIABLES DE ESTADO
c       svars(1:nfi) Deformaciones totales
c       svars(1+nfi:2*nfi) Esfuerzos generalizados
c       svars(1+2*nfi:2*nfi+nv) Variables internas con funcin de fluencia
c       svars(1+2*nfi+nv:2*nfi+nv+na) Variables internas sin funcin de fluencia
c       svars(1+2*nfi+nv+na:2*nfi+nv+na+ny) Fuerzas termodinﬂmicas
c
c******************************************************************
c
        include 'aba_param.inc'
c
        Real(kind=8):: rhs(mlvarx,*),amatrx(ndofel,ndofel),props(*),
     1          svars(*),energy(8),coords(mcrd,nnode),u(ndofel),
     2          du(mlvarx,*),v(ndofel),a(ndofel),time(2),params(*),
     3          adlmag(mdload,*),ddlmag(mdload,*),
     4          predef(2,npredf,nnode)
c
         Real(kind=8):: dtime
         Real(kind=8):: pnewdt
         integer lflags(*),jprops(*),jdltyp(mdload,*)
         integer ndofel,nrhs,nsvars,nprops,mcrd,nnode,jtype,kstep,kinc
         integer jelem,mdload,npredef,mlvarx,njprop
c
c******************************************************************
c                DECLARACION DE VARIABLES MDC
c******************************************************************
c
c        MATRIZ DE PROPIEDADES DEL ELEMENTO
c        prop(1)= longitud del elemento
c        prop(22)= masa total del elemento
c        prop(2)= seno del angulo que forma el elemento
c        prop(3)= coseno del angulo que forma el elemento
c        prop(4:n_prop)= propiedades restantes del elemento
c
c        desp= desplazamientos nodales al final del paso
c        velo= velocidades nodales al final del paso
c        acel= aceleraciones nodales al final del paso
c        param= parﬂmetros de la integracin numârica
c        opcion= opciones del calculo (dinﬂmico o estﬂtico,
c                        peque±os o grandes desplazamientos, etc.)
c        dtiempo= incremento de tiempo en el paso
c        varint_v_cero= var.int. con funcin de fluencia al principio del paso
c        varint_v = var.int.con funcin de fluencia al final del paso
c        varint_a_cero= var.int. sin funcin de fluencia al principio del paso
c        varint_a = var.int. sin funcin de fluencia al final del paso
c
c        residu= vector residual
c        jacob= jacobiano global
c        coord= coordenadas de los nudos del elemento
c        esf=esfuerzos al final del paso
c        esf_cero= esfuerzos al principio del paso
c        deftot=deformaciones totales al final del paso
c        deftot_cero= deformaciones totales al principio del paso
c        fterm= Fuerzas termodinamicas asociadas a las variables
c                     internas al final del paso
c        fterm_cero= Fuerzas termodinamicas al principio del paso
c        nq= nîmero de grados de libertad
c        nc= nîmero de coordenadas que definen cada nudo
c        nnu=nîmero de nudos del elemento
c        nv= nîmero de variables internas con funcin de fluencia
c        na= nîmero de variables internas sin funcin de fluencia
c        ny= nîmero de fuerzas termodinﬂmicas
c        n_prop= nîmero de propiedades del elemento
c        nfi= dimensin de la matriz de deformaciones y de esfuerzos
c        max_prop= tama±o de la matriz de propiedades
c        max_var_v= tama±o de la matriz de var.int. con funcin de fluencia
c        max_var_a= tama±o de la matriz de var.int. sin funcin de fluencia
c        max_opcion= tama±o de la matriz de opciones
c        max_param= tama±o de la matriz de parﬂmetros
c        max_desp= numero mﬂximo de desplazamientos
c        max_nud= mﬂximo nîmero de nudos del elemento
c        max_coo= mﬂximo nîmero de coordenadas por nodo
c        max_def= tama±o de la matriz de deformaciones
c        max_fterm= mﬂximo nîmero de fuerzas trmodinﬂmicas
c
c*********************************************************************
c
        integer max_prop,max_jprop,max_var_v,max_var_a
        integer max_desp,max_opcion,max_param
        integer max_nud,max_coo,max_def,max_fterm
        integer max_act
c
        parameter(max_prop=61)
        parameter(max_jprop=12)
        parameter(max_var_v=6)
        parameter(max_var_a=2)
        parameter(max_desp=6)
        parameter(max_opcion=3)
        parameter(max_param=17)
c
        parameter(max_nud=2)
        parameter(max_coo=2)
        parameter(max_def=3)
        parameter(max_fterm=4)
c
        parameter(max_act=6)
        Real(kind=8):: prop(max_prop),desp(max_desp),
     1             velo(max_desp),acel(max_desp),param(max_param),
     2             deftot_cero(max_def),esf_cero(max_def),
     2             dtiempo,deftot(max_def),esf(max_def)
c
        character(len=8):: opcion(max_opcion)
c
        integer n_prop,n_jprop,nq,jprop(max_jprop)
c
        Real(kind=8):: varint_v(max_var_v),varint_a(max_var_a),
     2         varint_v_cero(max_var_v),varint_a_cero(max_var_a),
     1         jacob(max_desp,max_desp),residu(max_desp),
     3         coord(max_coo,max_nud),tiempo_total
c
        Real(kind=8):: fterm(max_fterm),fterm_cero(max_fterm)
c
        Real(kind=8):: masa(max_desp,max_desp)
c
        logical activas(max_act),convergencia
	  character*6 damagefile 
c
        Real(kind=8):: t_ini_i,phi_i,io_i,iod_I_i,di_i,Rel_i,
     1                 t_ini_j,phi_j,io_j,iod_I_j,di_j,Rel_j,
     2                 corr_cero(2),corr(2),fcc_cero(2),fccorr(2),aux,
     3                 DhDq(2)
c
c*****************************************************************
c                TRADUCCION ABAQUS-MDC
c*****************************************************************
c
        tiempo_total = time(2)
c
	write (damagefile,'(I0)') jelem
	open(unit=100,file='C:\temp\DAMAGE'//trim(damagefile),
     1     access='append')
	write(100,20)kstep,time(2),svars(9),svars(10),svars(11),svars(12)
   20   format (I3,5(E12.5))
   	close(unit=100)
c	
        call trad_abaqus_mdc(props,jprops,coords,u,v,a,params,
     e                  jtype,lflags,svars,
     e                  dtime,nsvars,nprops,njprop,
     e                  mcrd,nnode,ndofel,pnewdt,
     e                  max_prop,max_jprop,max_var_v,max_var_a,max_desp,
     e                  max_opcion,max_param,max_fterm,
     e                  max_nud,max_coo,max_def,
     s                  coord,corr_cero,fcc_cero,
     s                  esf_cero,deftot_cero,fterm_cero,
     s                  prop,jprop,desp,velo,acel,param,opcion,
     s                  varint_v_cero,varint_a_cero,
     s                  dtiempo,nq,n_prop,n_jprop)
c
        if(opcion(1).eq."static".or.opcion(1).eq."dynamic")then
c
c******************************************************************
c        LLAMADA A MDC1
c        CALCULO DE LAS FUERZAS RESIDUALES, VARIABLES INTERNAS
c        Y JACOBIANO EN COORDENADAS GLOBALES AL
c        FINAL DEL PASO DE INTEGRACION
c        
c        ELEMENTO ELASTOPLASTICO DEGRADABLE 2D        
c******************************************************************
c
          write(*,*)jelem,kstep
          call mdc(prop,jprop,desp,velo,acel,param,opcion,
     e                  dtiempo,nq,n_jprop,n_prop,max_act,
     e                  max_prop,max_var_v,max_var_a,max_desp,
     e                  max_opcion,max_param,max_jprop,
     e                  varint_v_cero,varint_a_cero,
     e                  max_nud,max_coo,max_fterm,coord,
     e                  max_def,deftot_cero,esf_cero,
     e                  corr_cero,fcc_cero,
     s                  varint_v,varint_a,deftot,esf,
     s                  fterm,masa,DhDq,
     s                  residu,jacob)
c
c**********************************************************
c      C·lculo del grado de Corrosi¢n en las r¢tulas i y j
c**********************************************************
c
c      
c            Extremo en la RÛtula i
c
      t_ini_i=prop(56)
	phi_i=prop(57)
	io_i=prop(58)
	iod_I_i=prop(59)
	di_i=prop(60)
	Rel_i=prop(61) 	
	aux = (varint_v(3)+varint_v(5))/2.
c
      call cal_cor(
     e         varint_v(1),aux,dtiempo,tiempo_total,t_ini_i,phi_i,io_i,
     e         iod_I_i,di_i,Rel_i,corr_cero(1),DhDq(1),
     s         corr(1),fccorr(1))
c
c
c      
c            Extremo en la RÛtula j
c
      t_ini_j=prop(56)
	phi_j=prop(57)
	io_j=prop(58)
	iod_I_j=prop(59)
	di_j=prop(60)
	Rel_j=prop(61)
      aux = (varint_v(4)+varint_v(6))/2.
c
      call cal_cor(
     e         varint_v(2),aux,dtiempo,tiempo_total,t_ini_j,phi_j,io_j,
     e         iod_I_j,di_j,Rel_j,corr_cero(2),DhDq(2),
     s         corr(2),fccorr(2))
c
c
c ********************************************************************
c
        else if(opcion(1).eq."vel_jump")then
c
                call cal_def1(prop,max_def,max_desp,max_opcion,opcion,
     e                      desp,max_prop,
     e                      max_nud,max_coo,coord,deftot_cero,
     s                      deftot)
c
                call cal_masa(prop,max_prop,max_coo,max_nud,
     e                  max_desp,max_def,coord,
     s                  masa)
c
        else if(opcion(1).eq."ini_acel")then
c
                call cal_def1(prop,max_def,max_desp,max_opcion,opcion,
     e                      desp,max_prop,
     e                      max_nud,max_coo,coord,deftot_cero,
     s                      deftot)
c
                call cal_prop(prop,max_prop,deftot,max_def,esf_cero,
     e                    corr_cero,max_jprop,jprop,param,max_param,
     s                    DhDq)
c
                call esfuerzos(max_prop,max_jprop,prop,jprop,max_opcion,
     e              opcion,n_jprop,n_prop,max_var_v,max_var_a,max_param,
     e              param,varint_v_cero,varint_a_cero,max_def,max_fterm,
     e              deftot_cero,deftot,esf_cero,max_act,dtiempo,
     s              varint_v,varint_a,esf,fterm,activas,
     s                      convergencia)
c
                call cal_residu(max_prop,prop,max_desp,desp,max_opcion,
     e                  opcion,param,max_param,esf_cero,
     e                  max_def,esf,acel,velo,
     e                  nq,coord,max_nud,max_coo,
     s                  residu)
c
                call cal_masa(prop,max_prop,max_coo,max_nud,
     e                  max_desp,max_def,coord,
     s                  masa)
c
        end if
c
c******************************************************************
c       TRADUCCION MDC-ABAQUS
c******************************************************************
c
        call trad_mdc_abaqus(varint_v,varint_a,residu,jacob,
     e                     ndofel,nrhs,mlvarx,
     e                     max_var_v,max_var_a,max_desp,max_def,
     e                     max_fterm,max_jprop,esf,deftot,
     e                     max_param,param,fterm,jprop,
     e                     opcion,max_opcion,masa,corr,fccorr,
     s                     pnewdt,svars,rhs,amatrx)
c
c******************************************************************
c
        return
        end
c
c*****************************************************************
c
c        Elemento finito con variables internas para MDC
c        
c        MDC
c        VERSION 1.0
c
c******************************************************************
c
c        SUBRUTINA  DE  TRADUCCION  ABAQUS  --->  MDC
c        
c        ESCRITA POR Marça Eugenia Marante y Julio Flrez Lpez
c        Diembre 2001
c
c******************************************************************
c
        subroutine trad_abaqus_mdc(props,jprops,coords,u,v,a,params,
     e                  jtype,lflags,svars,
     e                  dtime,nsvars,nprops,njprop,
     e                  mcrd,nnode,ndofel,pnewdt,
     e                  max_prop,max_jprop,max_var_v,max_var_a,max_desp,
     e                  max_opcion,max_param,max_fterm,
     e                  max_nud,max_coo,max_def,
     s                  coord,corr_cero,fcc_cero,
     s                  esf_cero,deftot_cero,fterm_cero,
     s                  prop,jprop,desp,velo,acel,param,opcion,
     s                  varint_v_cero,varint_a_cero,
     s                  dtiempo,nq,n_prop,n_jprop)
c
       IMPLICIT NONE
c*************************************************************
c                DECLARACION DE LAS VARIABLES ABAQUS
c************************************************************
c
        integer nprops,njprop,nsvars,mcrd,nnode,ndofel
c
        integer jtype,lflags(4),jprops(*)
c
        Real(kind=8):: svars(*),props(*),
     2          coords(mcrd,nnode),u(ndofel),v(ndofel),
     3          a(ndofel),params(*)
     
               Real(kind=8):: dtime,pnewdt
c
c
c**********************************************************************
c                DECLARACION DE VARIABLES MDC
c**********************************************************************
c
      integer max_prop,max_jprop,max_var_v,max_var_a,max_desp,max_opcion
c
      integer max_param,max_nud,max_coo,max_def,max_fterm
c
      Real(kind=8):: prop(max_prop),desp(max_desp),
     1           velo(max_desp),acel(max_desp),param(max_param),
     2           deftot_cero(max_def),esf_cero(max_def)
c
      character(len=8):: opcion(max_opcion)
c
      integer nq,n_prop,n_jprop,jprop(max_jprop)
c
      Real(kind=8):: varint_v_cero(max_var_v),varint_a_cero(max_var_a),
     1       coord(max_coo,max_nud),fterm_cero(max_fterm),dtiempo,
     2       corr_cero(2),fcc_cero(2)
c
c
c**********************************************************************
c                DECLARACIONES LOCALES
c**********************************************************************
c
        integer nc,nfi,nv,na,ny,nnu,nt,i,j
c       Real(kind=8):: Es,kcor,aux
c**********************************************************************
c           TRADUCCION DE LAS COORDENADAS DEL ELEMENTO
c**********************************************************************
c
        n_jprop=njprop
c
        do i=1,n_jprop
            jprop(i)=jprops(i)
        end do
c
        nc=jprop(3)
        nnu=nnode
c
        do i=1,nc
                do j=1,nnu
                     coord(i,j)=coords(i,j)
                end do
        end do
c
c***********************************************************************
c        TRADUCCION DE LAS PARAMETROS DE IDENTIFICACION DEL ELEMENTO
c***********************************************************************
c
        nq=ndofel
        nfi=jprop(4)
        nv=jprop(5)
        na=jprop(6)
        ny=jprop(7)
        nt=nfi+nv+na
c
        n_prop=nprops-na-nv
c
        if(max_prop.lt.n_prop)then
             stop "error en la dimension del arreglo props"
        end if
c
c***********************************************************************
c        TRADUCCION DE LOS PARAMETROS DE INTEGRACION NıMERICA
c**********************************************************************
c
        param(1)=params(1)
        param(2)=params(2)
        param(3)=params(3)
        param(4)=jprops(1)
        param(5)=pnewdt
        param(6)=jprops(12)
c
c======================================================================
c
        do i=1,nfi+nv-4
                param(6+i)=props(i)
                if(param(6+i).le.0.)then
                    stop
     1                "error en los ceros numericos de la ley de estado
     2                  o de las funciones de plasticidad"
                end if
        end do
c
        do i=nfi+nv-3,nfi+nv
                param(6+i)=props(i)
                if(param(6+i).lt.0.)then
                    stop
     1                "error en los ceros numericos de las funciones
     2                  de danio"
                end if
        end do
c
        do i=nfi+nv+1,nfi+nv+na
                param(6+i)=props(i)
                if(param(6+i).le.0.)then
                    stop
     1                "error en los ceros numericos de las funciones
     2                  de plasticidad maxima"
                end if
        end do
c        
c=======================================================================
c
c***********************************************************************
c        TRADUCCION DE LAS PROPIEDADES DEL ELEMENTO
c***********************************************************************
c
        do i=1,n_prop-nfi
            prop(i+3)=props(i+nt)
        end do
c
c***********************************************************************
c        TRADUCCION DE LOS DESPLAZAMIENTOS NODALES, VELOCIDADES Y 
c        ACELERACIONES
c***********************************************************************
c
        do i=1,nq
                desp(i)=u(i)
                velo(i)=v(i)
                acel(i)=a(i)
        end do
c
c************************************************************************
c        TRADUCCION DE LAS OPCIONES
c  
c************************************************************************
c
        if (lflags(1).eq.1.or.lflags(1).eq.2)then
                opcion(1)="static"
        else if((lflags(1).eq.11.or.lflags(1).eq.12).and.
     1          (lflags(3).eq.1.or.lflags(3).eq.5))then
                opcion(1)="dynamic"
        else if(lflags(1).eq.11.or.lflags(1).eq.12.and.
     1          lflags(3).eq.4)then
                opcion(1)="vel_jump"
c        
        else if(lflags(1).eq.11.or.lflags(1).eq.12.and.
     1          lflags(3).eq.6)then
                opcion(1)="ini_acel"
        else
                stop "no es ni calculo estatico ni dinamico"
        end if
c
        if (lflags(2).eq.0)then
                opcion(2)="small"
        else if (lflags(2).eq.1)then
                opcion(2)="large"
        else
                stop "no son ni peque~as ni grandes deformaciones"
        end if
c
c************************************************************************
c        TRADUCION DE LAS VARIABLES INTERNAS
c************************************************************************
c
        do i=1,nfi
                deftot_cero(i)=svars(i)
                esf_cero(i)=svars(i+nfi)
        end do
c
        do i=1,nv
                varint_v_cero(i)=svars(i+2*nfi)
        end do
c        
        do i=1,na
                varint_a_cero(i)=svars(i+2*nfi+nv)
        end do
c
        do i=1,ny
                fterm_cero(i)=svars(i+2*nfi+nv+na)
        end do
c        
        corr_cero(1)=svars(ny+2*nfi+nv+na+1)
        corr_cero(2)=svars(ny+2*nfi+nv+na+2)
        fcc_cero(1)=svars(ny+2*nfi+nv+na+3)
        fcc_cero(2)=svars(ny+2*nfi+nv+na+4)
c
c*******************************************************************************************************************
c
c***********************************************************************
c        TRADUCCION DEL INCREMENTO DE TIEMPO
c***********************************************************************
c
        if(lflags(3).eq.1)then
                dtiempo=dtime
        else if(lflags(3).eq.5)then
                dtiempo=dtime/2.
        end if
c
c***********************************************************************
c
        return
        end
c
c*****************************************************************
c
c        Elemento finito con variables internas para MDC
c        
c        MDC
c        VERSION 1.0
c
c******************************************************************
c
c        SUBRUTINA  DE  TRADUCCION  MDC  --->  ABAQUS
c        
c        ESCRITA POR Marça Eugenia Marante y Julio Flrez Lpez
c        Febrero 2001
c
c******************************************************************
c
        subroutine trad_mdc_abaqus(varint_v,varint_a,residu,jacob,
     e                          ndofel,nrhs,mlvarx,
     e                          max_var_v,max_var_a,max_desp,max_def,
     e                          max_fterm,max_jprop,esf,deftot,
     e                          max_param,param,fterm,jprop,
     e                          opcion,max_opcion,masa,corr,fccorr,
     s                          pnewdt,svars,rhs,amatrx)
c
       IMPLICIT NONE
c*************************************************************
c                DECLARACION DE LAS VARIABLES ABAQUS
c************************************************************
c
        integer ndofel,nrhs,mlvarx
c
        Real(kind=8):: svars(1),rhs(mlvarx,*),amatrx(ndofel,ndofel),
     e               pnewdt
c
c
c******************************************************************
c                DECLARACION DE VARIABLES MDC 
c******************************************************************
c
        integer max_var_v,max_var_a,max_desp,max_def,max_opcion
c
        integer max_jprop,max_fterm,max_param
c
        Real(kind=8):: varint_v(max_var_v), varint_a(max_var_a)
        Real(kind=8):: jacob(max_desp,max_desp),param(max_param)
        Real(kind=8):: esf(max_def),deftot(max_def),residu(max_desp)
c
        integer jprop(max_jprop)
c
        Real(kind=8):: fterm(max_fterm),corr(2),fccorr(2)
c
        character(len=8):: opcion(max_opcion)
c
        Real(kind=8):: masa(max_desp,max_desp)
c
c******************************************************************
c                DECLARACION DE VARIABLES LOCALES
c******************************************************************
c
        integer nq,nfi,nv,na,ny,i,j
c
c**********************************************************************
c
        nq=ndofel
        nfi=jprop(4)
        nv=jprop(5)
        na=jprop(6)
        ny=jprop(7)
c
        pnewdt=param(5)
c
c**************************************************
c
        if(opcion(1).eq."static".or.opcion(1).eq."dynamic")then
           do i=1,nq
                rhs(i,1)=-residu(i)
                do j=1,nq
                         amatrx(i,j)=jacob(i,j)
                end do
           end do
c
           do i=1,nfi
                svars(i)=deftot(i)
                svars(i+nfi)=esf(i)
           end do
c
           do i=1,nv
                svars(i+2*nfi)=varint_v(i)
           end do
c
           do i=1,na
                svars(i+2*nfi+nv)=varint_a(i)
           end do 
c    
           do i=1,ny
                svars(i+2*nfi+nv+na)=fterm(i)
           end do         
c
           svars(ny+2*nfi+nv+na+1)=corr(1)
           svars(ny+2*nfi+nv+na+2)=corr(2)
           svars(ny+2*nfi+nv+na+3)=fccorr(1)
           svars(ny+2*nfi+nv+na+4)=fccorr(2)
c
c**************************************************
c
       else if(opcion(1).eq."vel_jump")then
           do i=1,nq
                do j=1,nq
                    amatrx(i,j)=masa(i,j)
                end do
           end do
c
c**************************************************
c
       else if(opcion(1).eq."ini_acel")then
          do i=1,nq
                  rhs(i,1)=-residu(i)
                do j=1,nq
                   amatrx(i,j)=masa(i,j)
                end do
          end do
c
c**************************************************
c
       else
            stop "error en opcion(1)"
c
       end if
c
         pnewdt=param(5)
c
       return
       end
c
c******************************************************************
c
c        SUBRUTINA  MDC1
c        CALCULO DE LAS DEFORMACIONES TOTALES, ESFUERZOS,
c        VARIABLES INTERNAS, FUERZAS  TERMODINAMICAS,
c        JACOBIANO GLOBAL Y
c        FUERZAS RESIDUALES AL FINAL DEL PASO
c        
c        VERSION 1.0
c        ESCRITA POR Marça Eugenia Marante y Julio Flrez Lpez
c        Marzo 2001
c
c******************************************************************
c
        subroutine mdc(prop,jprop,desp,velo,acel,param,opcion,
     e                        dtiempo,nq,n_jprop,n_prop,max_act,
     e                        max_prop,max_var_v,max_var_a,max_desp,
     e                        max_opcion,max_param,max_jprop,
     e                        varint_v_cero,varint_a_cero,
     e                        max_nud,max_coo,max_fterm,coord,
     e                        max_def,deftot_cero,esf_cero,
     e                        corr_cero,fcc_cero,
     s                        varint_v,varint_a,deftot,esf,
     s                        fterm,masa,DhDq,
     s                        residu,jacob)
c
       IMPLICIT NONE
c******************************************************************
c                DECLARACION DE VARIABLES MDC
c******************************************************************
c
c        PROP: MATRIZ DE PROPIEDADES DEL ELEMENTO 
c        prop(1)= longitud del elemento
c        prop(2:n_prop)= propiedades restantes del elemento
c
c     PAR¬METROS DE IDENTIFICACI‘N DEL ELEMENTO
c     jprop(1)=nc:nîmero de coordenadas que definen cada nudo
c     jprop(2)=nfi:dimensin de la matriz de deformaciones y de esfuerzos
c     jprop(3)=nv:nîmero de variables internas con funcin de fluencia
c     jprop(4)=na:nîmero de variables internas sin funcin de fluencia
c     jprop(5)=ny:nîmero de fuerzas termodinﬂmicas
c
c        PAR¬METROS DE INTEGRACI‘N NUM–RICA
c        param(1:3) Coeficientes del mâtodo de newmark
c        param(4)=max_ite=maximo No.de iteraciones
c        param(5)= max_paso=maximo No.de pasos de integracion 
c        param(7:6+nfi+nv+na)= Ceros numâricos de {R},{T},{G}
c
c        desp= desplazamientos nodales al final del paso
c        velo= velocidades nodales al final del paso
c        acel= aceleraciones nodales al final del paso
c        param= parﬂmetros de la integracin numârica
c        opcion= opciones del calculo (dinﬂmico o estﬂtico,
c                        peque±os o grandes desplazamientos, etc.)
c        dtiempo= incremento de tiempo en el paso
c        varint_v_cero= var.int. con funcin de fluencia al principio del paso
c        varint_v = var.int.con funcin de fluencia al final del paso
c        varint_a_cero= var.int. sin funcin de fluencia al principio del paso
c        varint_a = var.int. sin funcin de fluencia al final del paso
c        residu= vector residual
c        jacob= jacobiano global
c        coord= coordenadas de los nudos del elemento
c        esf=esfuerzos al final del paso
c        esf_cero= esfuerzos al principio del paso
c        deftot=deformaciones totales al final del paso
c        deftot_cero= deformaciones totales al principio del paso
c        fterm= Fuerzas termodinamicas asociadas a las variables 
c                     internas al final del paso
c        fterm_cero= Fuerzas termodinamicas al principio del paso
c        nq= nîmero de grados de libertad
c        nnu=nîmero de nudos del elemento
c        n_prop= nîmero de propiedades del elemento
c     n_jprop=nîmero de parﬂmetros de identificacin del elemento
c        
c        max_prop= tama±o de la matriz de propiedades
c     max_jprop= tama±o de parﬂmetros de identificacin del elemento
c        max_var_v= tama±o de la matriz de var.int. con funcin de fluencia
c        max_var_a= tama±o de la matriz de var.int. sin funcin de fluencia 
c        max_opcion= tama±o de la matriz de opciones
c        max_param= tama±o de la matriz de parﬂmetros
c        max_desp= numero mﬂximo de desplazamientos
c        max_nud= mﬂximo nîmero de nudos del elemento
c        max_coo= mﬂximo nîmero de coordenadas por nodo=2
c        max_def= tama±o de la matriz de deformaciones
c        max_fterm= mﬂximo nîmero de fuerzas termodinﬂmicas
c
c******************************************************************
c
        integer max_prop,max_jprop,max_var_v,max_var_a,max_desp
c
        integer max_opcion,max_param,max_nud,max_coo,max_def,max_fterm
c
        integer max_act
c
        Real(kind=8):: prop(max_prop),desp(max_desp),
     1          velo(max_desp),acel(max_desp),param(max_param),dtiempo
c
        character(len=8):: opcion(max_opcion)
c
        integer nq,n_prop,n_jprop,jprop(max_jprop)
c
        Real(kind=8):: varint_v(max_var_v),varint_a(max_var_a),
     1             varint_v_cero(max_var_v),varint_a_cero(max_var_a),
     2             jacob(max_desp,max_desp),
     3             coord(max_coo,max_nud),masa(max_desp,max_desp)
c
        Real(kind=8):: deftot_cero(max_def),esf_cero(max_def),
     1         deftot(max_def),esf(max_def),residu(max_desp)
c
        Real(kind=8):: fterm(max_fterm),corr_cero(2),fcc_cero(2),
     1                 DhDq(2)
c
c******************************************************************
c                DECLARACIONES LOCALES
c******************************************************************
c
c        activas= matriz que indica el estado de todas las variables
c                internas 
c        max_act= tama~o maximo de la matriz activas
c        deftot_int= matriz de deformaciones totales intermedias.
c                                Se usa para la gestion de paso local
c        max_deftot_int= tama±o maximo de deftot_int.
c                                        Debe ser igual a max_def
c        alpha= parametro que indica el numero de subintervalos en los 
c                        que se divide cada paso
c
c******************************************************************
c
       integer max_deftot_int,i
       parameter(max_deftot_int=3)
       logical activas(max_act),convergencia
       Real(kind=8):: alpha,deftot_int(max_deftot_int),
     1                prop_II(max_prop)
c
c******************************************************************
c        CALCULO DEL VECTOR DE DEFORMACIONES TOTALES
c        AL FINAL DEL INCREMENTO
c******************************************************************
c
       call cal_def1(prop,max_def,max_desp,max_opcion,opcion,
     e                desp,max_prop,
     e                max_nud,max_coo,coord,deftot_cero,
     s                deftot)
c
c*********************************************************************
c        Lazo para la gestion del paso local
c*********************************************************************
c
       alpha=1.0        
       convergencia=.false.
       do while (.not.convergencia.or.alpha.ne.1.)
c                
c*********************************************************************
c        calculo de las deformaciones intermedias para la gestion
c        del paso local
c*********************************************************************
c
          do i=1,max_def
             deftot_int(i)=alpha*deftot(i)+(1.-alpha)*deftot_cero(i)
          end do
c
c******************************************************************
c        CALCULO DE LOS COEFICIENTES DEL MODELO
c******************************************************************
c
          do i=1,max_prop
		    prop_II(i)=prop(i)
	    end do
c
          call cal_prop(prop,max_prop,deftot,max_def,esf_cero,
     e              corr_cero,max_jprop,jprop,param,max_param,
     s              DhDq)
c
c******************************************************************
c        CALCULO DE LOS ESFUERZOS, LAS VARIABLES INTERNAS
c        Y LAS FUERZAS TERMODINAMICAS AL FINAL DEL PASO
c******************************************************************
c
          call esfuerzos(max_prop,max_jprop,prop,jprop,max_opcion,
     e            opcion,n_jprop,n_prop,max_var_v,max_var_a,max_param,
     e            param,varint_v_cero,varint_a_cero,max_def,max_fterm,
     e            deftot_cero,deftot_int,esf_cero,max_act,dtiempo,
     s            varint_v,varint_a,esf,fterm,activas,
     s                    convergencia)
c
c******************************************************************
c         GESTION DEL PASO LOCAL
c******************************************************************
c
          call cambio_de_paso(max_def,max_deftot_int,max_var_v,
     e        max_var_a,deftot,deftot_int,varint_v,varint_a,esf,
     s        varint_v_cero,varint_a_cero,esf_cero,max_param,
     s              param,convergencia,deftot_cero,alpha)
c
          if (.not.convergencia) then
		   do i=1,max_prop
		        prop(i)=prop_II(i)
	       end do
	    end if
c
       end do
c
c******************************************************************
c        CALCULO DE LAS FUERZAS RESIDUALES INTERNAS
c******************************************************************
c
       call cal_residu(max_prop,prop,max_desp,desp,max_opcion,
     e                  opcion,param,max_param,esf_cero,
     e                  max_def,esf,acel,velo,
     e                  nq,coord,max_nud,max_coo,
     s                  residu)
c
c******************************************************************
c        CALCULO DEL JACOBIANO GLOBAL
c******************************************************************
c
       call cal_jacob(max_prop,prop,max_desp,desp,max_opcion,
     e           max_jprop,jprop,opcion,nq,n_prop,max_var_v,
     e           max_var_a,max_def,max_coo,max_nud,esf,dtiempo,
     e           coord,deftot,varint_v,varint_a,
     e           param,max_param,max_act,activas,
     e           varint_v_cero,varint_a_cero,deftot_cero,esf_cero,
     s           jacob)
c
c******************************************************************
c        CALCULO DE LA MATRIZ DE MASA
c******************************************************************
c
       call cal_masa(prop,max_prop,max_coo,max_nud,
     e                  max_desp,max_def,coord,
     s                  masa)
c
c******************************************************************
c                
        return
        end
c
c******************************************************************
c        SUBRUTINA  CAL_DEF1
c        CALCULO DE LAS DEFORMACIONES TOTALES

c
c        VERSION 1.0
c        ESCRITA POR Marça Eugenia Marante y Julio Flrez Lpez
c        Modificada por Nayive Jaramillo 
c        (Enero 2002)
c
c******************************************************************
c
        subroutine cal_def1(prop,max_def,max_desp,max_opcion,opcion,
     e                      desp,max_prop,
     e                      max_nud,max_coo,coord,deftot_cero,
     s                      deftot)
c
       IMPLICIT NONE
c******************************************************************
c        DECLARACION DE VARIABLES MDC
c******************************************************************
c
c        prop= matriz de propiedades del elemento
c                  en peque~nas deformaciones        
c        prop(1)= longitud del elemento
c        prop(2)=seno del angulo con la horizontal
c       prop(3)=coseno del angulo con la horizontal
c                        en Calculos Dinamicos
c       prop(34)= masa total del elemento
c       prop(4:n_prop)= propiedades restantes del elemento
c
c        desp= desplazamientos nodales al final del paso
c        max_prop= tama±o de la matriz de propiedades
c        max_opcion= tama±o de la matriz de opciones
c        opcion=opciones del calculo (dinamico o estatico,
c       peque~os o grandes desplazamientos, etc)
c        opcion(1)= "static" o "dynamic"
c        opcion(2)= "small" o "large"
c        max_desp= numero maximo de desplazamientos
c        coord= coordenadas de los nudos del elemento
c        max_nud= maximo numero de nudos del elemento
c        max_coo= maximo numero de coordenadas
c        max_def= tama±o de la matriz de deformaciones
c        deftot=deformaciones totales al final del paso
c        deftot_cero= deformaciones totales al principio del paso
c
c******************************************************************
c
        integer max_prop,max_desp,max_opcion
c
        integer max_nud,max_coo,max_def
c
        Real(kind=8):: prop(max_prop),desp(max_desp)
c
        Real(kind=8):: coord(max_coo,max_nud)
c
        Real(kind=8):: deftot(max_def),deftot_cero(max_def)
c
        character(len=8):: opcion(max_opcion)
c
c*********************************************************************
c               DECLARACIONES LOCALES
c*********************************************************************
c
c       b= matriz de transformacion desp-deftot
c       long= longitud del elemento
c       seno= seno del angulo con la horizontal
c       coseno= coseno del angulo con la horizontal
c       b= matriz de transformacion desp-deftot
c       alfacero=angulo de inclinacion de la barra al principio del paso
c       alfa=angulo de inclinacion de la barra al final del paso
c       X()=coordenada horizontal al final del paso para grandes desp.  
c       Y()=coordenada vertical al final del paso para grandes desp.    
c
c*********************************************************************
c
        Real(kind=8):: long,seno,coseno,b(3,6),alfacero,alfa,longt
        Real(kind=8):: x1,x2,y1,y2
        Real(kind=8):: PI
c
c*********************************************************************
        if (opcion(2).eq."small")then
c
c               CALCULO EN PEQUE~AS DEFORMACIONES
c
c               calculo de la longitud del elemento
c
                prop(1)=dsqrt((coord(1,2)-coord(1,1))**2+
     1                  (coord(2,2)-coord(2,1))**2)
c
                long=prop(1)
c
c               CALCULO DEL SENO
C
                prop(2)=(coord(2,2)-coord(2,1))/prop(1)
                seno=prop(2)
C
C               CALCULO DEL COSENO
c
                prop(3)=(coord(1,2)-coord(1,1))/prop(1)
                coseno=prop(3)
c
c               CALCULO DE LA MATRIZ DE TRANSFORMACION EN 
C               LA CONFIGURACION ORIGINAL
C
                call cal_b(long,seno,coseno,
     e                          max_def,max_desp,
     s                          b)
c
c               CALCULO DE LAS DEFORMACIONES
C
                call pro_2mat(max_def,max_desp,max_desp,1,
     e          max_def,1,
     e          max_def,max_desp,1,
     e          b,desp,
     s          deftot)
c
        else if(opcion(2).eq."large")then
c 
c ------------------------------------------------------------
c               OPCION PARA GRANDES DEFORMACIONES
c               ESCRITA POR: A. LOPEZ INOJOSA
c               REALIZADA EL 2/2/94
c
c
c              CALCULO DE LA LONGITUD INICIAL:
c
                long=dsqrt((coord(1,2)-coord(1,1))**2+
     1                  (coord(2,2)-coord(2,1))**2)
c
c               COORDENADAS AL FINAL DEL PASO:
c
                 x1=coord(1,1)+desp(1)
                 x2=coord(1,2)+desp(4)
                 y1=coord(2,1)+desp(2)
                 y2=coord(2,2)+desp(5)
c
c
c               CALCULO DE LA LONGITUD EN LA CONF. DEFORMADA:
                longt=dsqrt((x2-x1)**2+(y2-y1)**2)
c
c*******************************************************************
c               Calculo del seno y coseno en Configuracion Inicial
c               (Calculados para condicion logica):
c*******************************************************************
c
                seno=(coord(2,2)-coord(2,1))/long
                coseno=(coord(1,2)-coord(1,1))/long
c
c********************************************************************           
c
                PI=3.14159265359
c
c*******************************************************************
c               CALCULO DEL ANGULO EN LA CONFIGURACION INICIAL
c 
c               Para 90 Grados:
                if (seno.eq.1.)then
                  alfacero=PI/2.
c
c               Para 270 Grados
                else if(seno.eq.-1.)then
                  alfacero=3.*PI/2.
c
c               Primer y cuarto Cuadrante:              
                else if (coseno.gt.0.)then
                   alfacero=atan
     1               ((coord(2,2)-coord(2,1))/(coord(1,2)-coord(1,1)))
c
c               Segundo y Tercer Cuadrante
                else
                   alfacero=(atan
     1               ((coord(2,2)-coord(2,1))/(coord(1,2)-coord(1,1))))+
     2               PI
c
                end if              
c       
c
c***************************************************************************
c               CALCULO DE LAS PROPIEDADES DE LA BARRA PARA GRANDES DESP.
c
                   prop(1)=longt
c
c                 seno:
                   prop(2)=(y2-y1)/prop(1)
c                 coseno:
                   prop(3)=(x2-x1)/prop(1)
c
c***************************************************************************

c               CALCULO DEL ANGULO EN LA CONFIGURACION DEFORMADA
c 
c               Para 90 Grados:
                if (prop(2).eq.1.)then
                        alfa=PI/2.
c
c               Para 270 Grados
                else if(prop(2).eq.-1.)then
                        alfa=3.*PI/2.
c
c               Primer y Cuarto Cuadrante:              
                else if (prop(3).ge.0.)then
                        alfa=atan((y2-y1)/(x2-x1))
c
c               Segundo y Tercer cuadrante
                else
                        alfa=(atan((y2-y1)/(x2-x1)))+PI
c
                end if
c
c
c               CALCULO DE LAS DEFORMACIONES ANGULARES
c
                deftot(1)=desp(3)-(alfacero-alfa)
                deftot(2)=desp(6)-(alfacero-alfa)
c
c               CALCULO DE LA DEFORMACION AXIAL:
c
                deftot(3)=longt-long
c                               
c
c---------------------------------------------------------------------- 
        else
                stop
     1          "error en la matriz opcion"
        end if
c
        return
        end
c **************************************************************
c **************************************************************
c
c               SUBRUTINA CAL_B
c               CALCULO DE LA MATRIZ DE TRANSFORMACION 
c               ESCRITA POR: J. FLOREZ LOPEZ
c               COMENZADA EL 8/10/93
c                          Modificada en enero (2002)
c
c****************************************************************
c
        subroutine cal_b(long,seno,coseno,
     e                  max_def,max_desp,
     s                  b)
c
       IMPLICIT NONE
c******************************************************************
c               SIGNIFICADO DE LAS VARIABLES
c*****************************************************************
c
c       long=longitud del elemento
c       seno=seno del angulo con la horizontal
c       coseno= coseno del angulo con la horizontal
c       b= matriz de transformacion desp-deftot
c
c******************************************************************
c               DECLARACIONES
c*****************************************************************
c
        integer max_def,max_desp
        Real(kind=8):: b(max_def,max_desp)
c
        Real(kind=8):: seno,coseno,long
c
c*****************************************************************
c
c
        b(1,1)=seno/long
        b(1,2)=-coseno/long
        b(1,3)=1.
        b(1,4)=-seno/long
        b(1,5)=coseno/long
        b(1,6)=0.
        b(2,1)=seno/long
        b(2,2)=-coseno/long
        b(2,3)=0.
        b(2,4)=-seno/long
        b(2,5)=coseno/long
        b(2,6)=1.
        b(3,1)=-coseno
        b(3,2)=-seno
        b(3,3)=0.
        b(3,4)=coseno
        b(3,5)=seno
        b(3,6)=0.
c
c
        return
        end
c
c *****************************************************************
c                SUBRUTINA PRO_2MAT
c                PRODUCTO DE DOS MATRICES A x B = C
c                DIMENSIONES DE A: MAX_FIL_A,MAX_COL_A
c                FILAS Y COLUMNAS DE A: FIL_C,DIM
c                DIMENSIONES DE B: MAX_FIL_B,MAX_COL_B
c                FILAS Y COLUMNAS DE B: DIM,COL_C
c                DIMENSIONES DE C: MAX_FIL_C,MAX_COL_C
c                ESCRITA POR: J. FLOREZ LOPEZ
c******************************************************************
c
        subroutine pro_2mat(max_fil_a,max_col_a,max_fil_b,max_col_b,
     1          max_fil_c,max_col_c,
     1          fil_c,dim,col_c,a,b,
     1          c)
c
       IMPLICIT NONE
c******************************************************************
c                DECLARACIONES
c******************************************************************
c
        integer max_fil_a,max_col_a,max_fil_b,max_col_b
        integer fil_c,dim,col_c,max_fil_c,max_col_c
        integer i,j,k
c
        Real(kind=8):: a(max_fil_a,max_col_a),b(max_fil_b,max_col_b),
     1         c(max_fil_c,max_col_c)
c
c******************************************************************
c
        do i=1,fil_c
                do j=1,col_c
                        c(i,j)=0.
                        do k=1,dim
                                c(i,j)=c(i,j)+a(i,k)*b(k,j)
                        end do
                end do
        end do
c
        return
        end
c
c *****************************************************************
c                SUBRUTINA PRO_3MAT
c                PRODUCTO DE TRES MATRICES A x B x C = D
c                DIMENSIONES DE A: MAX_FIL_A,MAX_COL_A
c                FILAS Y COLUMNAS DE A: FIL_D,DIM1
c                DIMENSIONES DE B: MAX_FIL_B,MAX_COL_B
c                FILAS Y COLUMNAS DE B: DIM1,DIM2
c                DIMENSIONES DE C: MAX_FIL_C,MAX_COL_C
c                FILAS Y COLUMNAS DE C:DIM2,COL_D
c                DIMENSIONES DE D: MAX_FIL_D,MAX_COL_D
c        
c                ESCRITA POR: J. FLOREZ LOPEZ
c******************************************************************
c
        subroutine pro_3mat(max_fil_a,max_col_a,max_fil_b,max_col_b,
     1          max_fil_c,max_col_c,max_fil_d,max_col_d,
     1          fil_d,dim1,dim2,col_d,a,b,c,
     1          d)
c
       IMPLICIT NONE
c******************************************************************
c                DECLARACIONES
c******************************************************************
c
        integer max_fil_a,max_col_a,max_fil_b,max_col_b
        integer max_fil_d,max_col_d
        integer fil_d,dim1,dim2,col_d,max_fil_c,max_col_c
        integer i,j,k,l
c
        Real(kind=8):: a(max_fil_a,max_col_a),b(max_fil_b,max_col_b),
     1         c(max_fil_c,max_col_c),d(max_fil_d,max_col_d)
c
c******************************************************************
c
        do i=1,fil_d
               do j=1,col_d
                        d(i,j)=0.
                        do k=1,dim1
                           do l=1,dim2
                             d(i,j)=d(i,j)+a(i,k)*b(k,l)*c(l,j)
                           end do
                        end do
               end do
        end do
c
        return
        end
c
c******************************************************************
c        SUBRUTINA CAL_Bt
c        CALCULO DE LA MATRIZ DE TRANSFORMACION 
c        ESCRITA POR: Marça Eugenia Marante y Julio Flrez Lpez
c        (Abril 2001)
c
c******************************************************************
c
        subroutine cal_bt(b,max_esf,max_lib,
     s                                  bt)
c
       IMPLICIT NONE
c******************************************************************
c        SIGNIFICADO DE LAS VARIABLES
c******************************************************************
c
c        long=longitud del elemento
c        bt= matriz de transformacion desp-deftot
c        max_lib= max_desp
c        max_esf= max_def
c
c******************************************************************
c                DECLARACIONES
c******************************************************************
c
        integer max_lib,max_esf,i,j
        Real(kind=8):: b(max_esf,max_lib)
        Real(kind=8):: bt(max_lib,max_esf)
c
c******************************************************************        
c        
        do i=1,max_esf
                do j=1, max_lib
                         bt ( j, i) = b (i, j)
                end do
        end do
c
        return
        end
c
c******************************************************************       ****************************
c
c        SUBRUTINA CAL_PROP                                                     Generador 1
c        CALCULO DE LOS COEFICIENTES DEL MODELO HISTERETICO  DE
c        PORTICO DEGRADABLE  A  PARTIR  DE LOS DIAGRAMAS DE INTERACCI‘N             y
c        DE LOS COEFICIENTES
c
c     Escrita por Marça Eugenia Marante y Julio Flrez Lpez                    Generador 2
c     (Mayo 2001) 
c        Modificada por Nayive Jaramillo (ENERO 2002)
c
c        Modificada por Ricardo PicÛn, Scarlet Montilla y Julio FlÛrez LÛpez
c
c******************************************************************       ****************************
c
        subroutine cal_prop(prop,max_prop,deftot,max_def,esf_cero,
     e                  corr_cero,max_jprop,jprop,param,max_param,
     s                  DhDq)
c
       IMPLICIT NONE
c******************************************************************
c               DECLARACION DE VARIABLES NUCLEO
c
c******************************************************************
c
c       prop=matriz de propiedades del elemento
c       max_prop=tama±o de la matriz de propiedades
c       jprop = parﬂmetros de identificacin del elemento
c       max_jprop=tama±o de la matriz de parﬂmetros de identificacin del elemento
c       deftot=deformaciones totales al final del paso
c       max_def=tamano de la matriz de deformaciones
c
c******************************************************************
c
      integer max_prop,max_jprop,max_def,max_param
      integer jprop(max_jprop)
      Real(kind=8):: prop(max_prop),deftot(max_def),param(max_param)
      Real(kind=8):: esf_cero(max_def),corr_cero(2),DhDq(2)
c
c******************************************************************
c               DECLARACIONES LOCALES
c******************************************************************
c
c       fa = deformacion axial de la cuerda
c       n  = fuerza axial
c       c  = coeficiente del criterio de plasticidad
c       my = coeficiente del criterio de plasticidad                                  
c       q  = coeficiente del criterio de plasticidad
c       gcr = coeficiente del criterio de plasticidad
c       Los coeficientes del criterio de plasticidad se calculan
c       tanto condiciones positivas como negativas en i y j 
c******************************************************************
c
      Real(kind=8):: fa,n,l,ei,ea,ci,cj,h,rec,d
      Real(kind=8):: gcr_pos_i,q_pos_i,my_pos_i
	Real(kind=8):: c_pos_i,dp_pos_i,du_pos_i
      Real(kind=8):: gcr_neg_i,q_neg_i,my_neg_i
      Real(kind=8):: c_neg_i,dp_neg_i,du_neg_i
	Real(kind=8):: gcr_pos_j,q_pos_j,my_pos_j
      Real(kind=8):: c_pos_j,dp_pos_j,du_pos_j
	Real(kind=8):: gcr_neg_j,q_neg_j,my_neg_j
	Real(kind=8):: c_neg_j,dp_neg_j,du_neg_j,
     1               c_i_plus,q_i_plus,c_j_plus,q_j_plus,
     2               c_i_cero,q_i_cero,c_j_cero,q_j_cero,
	3               prop_II(max_prop),masa,masa_total,
     4               alfa_pos_i,alfa_neg_i,
	5               beta_pos_i,beta_neg_i
c
      integer i
c
c******************************************************************
c       Propiedades del elemento
c******************************************************************
c
      do i=1,max_prop
	   prop_II(i)=prop(i)
	end do
c
      l=prop_II(1)
      h=prop_II(5)
	rec=prop_II(18)
	ea=prop_II(39)*prop_II(4)*prop_II(5)
	ei=prop_II(39)*prop_II(4)*prop_II(5)*prop_II(5)*prop_II(5)/12.
	masa=prop_II(38)
	alfa_pos_i=prop_II(55)
      alfa_neg_i=prop_II(55)
      beta_pos_i=1.-alfa_pos_i
	beta_neg_i=1.-alfa_neg_i
c              
c******************************************************************
c       Alargamiento de la cuerda
c******************************************************************
c
      fa=deftot(3)
c
c******************************************************************
c       Fuerza axial 
c******************************************************************
c
        n=ea*fa/l/1000.
        ci=corr_cero(1)
        cj=corr_cero(2)
	  d=h-rec
c
c******************************************************************
c       Calculo de los coeficientes del criterio de plasticidad
c        para el modelo histeretico ci+ -,qi+ -,myi+ -,gcri+ -
c                   cj+ -,qj+ -,myj+ -,gcrj+ -
c
c******************************************************************
c      Coeficiente en i
c
       call Prop_acero_i(
     e            prop_II,ci)
c
        Call cal_generador_1(
     e       max_prop,prop_II,n,ei,l,esf_cero,
     s       gcr_pos_i,q_pos_i,my_pos_i,c_pos_i,dp_pos_i,du_pos_i,
     s       gcr_neg_i,q_neg_i,my_neg_i,c_neg_i,dp_neg_i,du_neg_i,
     s       gcr_pos_j,q_pos_j,my_pos_j,c_pos_j,dp_pos_j,du_pos_j,
     s       gcr_neg_j,q_neg_j,my_neg_j,c_neg_j,dp_neg_j,du_neg_j)
c
c
c
c******************************************************************
c******************************************************************
c       Coeficientes positivos del criterio de plasticidad en la
c       matriz de propiedades del elemento EXTREMO i
c******************************************************************
        prop(6)=c_pos_i
        prop(7)=q_pos_i
        prop(8)=my_pos_i
        prop(9)=gcr_pos_i
c******************************************************************
c       Coeficientes negativos del criterio de plasticidad en la              
c       matriz de propiedades del elemento EXTREMO i
c******************************************************************
        prop(14)=c_neg_i
        prop(15)=q_neg_i
        prop(16)=my_neg_i
        prop(17)=gcr_neg_i
c******************************************************************
c       Coeficiente en j
c
       call Prop_acero_j(
     e            prop_II,cj)
c
        Call cal_generador_1(
     e       max_prop,prop_II,n,ei,l,esf_cero,
     s       gcr_pos_i,q_pos_i,my_pos_i,c_pos_i,dp_pos_i,du_pos_i,
     s       gcr_neg_i,q_neg_i,my_neg_i,c_neg_i,dp_neg_i,du_neg_i,
     s       gcr_pos_j,q_pos_j,my_pos_j,c_pos_j,dp_pos_j,du_pos_j,
     s       gcr_neg_j,q_neg_j,my_neg_j,c_neg_j,dp_neg_j,du_neg_j)
c
c
c******************************************************************
c       Coeficientes positivos del criterio de plasticidad en la
c       matriz de propiedades del elemento EXTREMO j
c******************************************************************
        prop(10)=c_pos_j
        prop(11)=q_pos_j
        prop(12)=my_pos_j
        prop(13)=gcr_pos_j
c******************************************************************
c       Coeficientes negativos del criterio de plasticidad enla              
c       matriz de propiedades del elemento EXTREMO j
c******************************************************************
        prop(18)=c_neg_j
        prop(19)=q_neg_j
        prop(20)=my_neg_j
        prop(21)=gcr_neg_j
c******************************************************************
c       Coeficientes positivos y negativos para considerar
c       el efecto de pinching en los extremos i y j
c******************************************************************
c
	  prop(4)=ei/1000.
        prop(5)=ea/1000.
        masa_total=masa/1000.
	  prop(22)=masa_total
	  prop(23)=alfa_pos_i
	  prop(24)=alfa_neg_i
	  prop(25)=d
c
c******************************************************************
c
c       Determinacion de @h/@c    y   @h/@c
c          para calcular      @h/@q
c
c*******************************************************************
c
        ci = corr_cero(1) + 0.001
        cj = corr_cero(2) + 0.001
c
c          Extremo i
c
       call Prop_acero_i(
     e            prop_II,ci)
c
        Call cal_generador_1(
     e       max_prop,prop_II,n,ei,l,esf_cero,
     s       gcr_pos_i,q_pos_i,my_pos_i,c_pos_i,dp_pos_i,du_pos_i,
     s       gcr_neg_i,q_neg_i,my_neg_i,c_neg_i,dp_neg_i,du_neg_i,
     s       gcr_pos_j,q_pos_j,my_pos_j,c_pos_j,dp_pos_j,du_pos_j,
     s       gcr_neg_j,q_neg_j,my_neg_j,c_neg_j,dp_neg_j,du_neg_j)
c
c         Extremo j
c
       call Prop_acero_j(
     e            prop_II,cj)
c
        Call cal_generador_1(
     e       max_prop,prop_II,n,ei,l,esf_cero,
     s       gcr_pos_i,q_pos_i,my_pos_i,c_pos_i,dp_pos_i,du_pos_i,
     s       gcr_neg_i,q_neg_i,my_neg_i,c_neg_i,dp_neg_i,du_neg_i,
     s       gcr_pos_j,q_pos_j,my_pos_j,c_pos_j,dp_pos_j,du_pos_j,
     s       gcr_neg_j,q_neg_j,my_neg_j,c_neg_j,dp_neg_j,du_neg_j)
c
       c_i_plus = (c_pos_i + c_neg_i)/2.
       c_i_cero = (prop(6) + prop(14))/2.
       c_j_plus = (c_pos_j + c_neg_j)/2.
       c_j_cero = (prop(10) + prop(18))/2.
c
c
       q_i_plus = (q_pos_i + q_neg_i)/2.
       q_i_cero = (prop(7) + prop(15))/2.
       q_j_plus = (q_pos_j + q_neg_j)/2.
       q_j_cero = (prop(11) + prop(19))/2.
c
c
      DhDq(1) = (c_i_plus - c_i_cero) / (q_i_plus - q_i_cero)
      DhDq(2) = (c_j_plus - c_j_cero) / (q_j_plus - q_j_cero)
c
c
c
        return
        end
c
*******************************************************************
c
c        SUBRUTINA ESFUERZOS
c
c        CALCULO DE LOS ESFUERZOS,
c        LAS VARIABLES INTERNAS Y SUS FUERZAS TERMODINAMICAS
c        MEDIANTE LA RESOLUCION DEL SIGUIENTE SISTEMA DE
c        ECUACIONES NO LINEALES:
c
c        
c        R(esfuerzos,varint_v,varint_a,deftot)=0
c        T= (Vr - Vr-1)j fj(esfuerzos,varint_v,varint_a,deftot)=0 donde f<=0
c        G(esfuerzos,varint_v,varint_a,deftot)=0
c
c        ....
c        DONDE
c
c        R representa a la ley de estado
c        T representa la ley de evolucin de las variables internas Vi
c        G representa la Ley de Danio que define las variables internas Ai 
c        f es la funcion de fluencia de la variable interna Vj
c
c        Modificada por Nayive Jaramillo (Enero 2002)
c
c******************************************************************
c
c        Modificada por Lorena Suarez
c        Junio 2003
c        La modificacion consistio en realizar un manejo de errores para los casos en que
c        la rutina colapsaba debido a errores numericos (basicamente valores imaginarios)
c        Esto fue arreglado facilmente obligando a la rutina a reducir el paso a nivel local
c        Se modificaron las rutinas de esfuerzos, cal_f,cal_dano, y todas aquella rutinas que hacen uso
c        de cal_f. Como yo no sabia en que puntos o llamadas a cal_f se colgaban los programas
c        mas aun teniendo en cuenta que pueden ser de diversos tipos, realice verificaciones de
c        los valores arrojados por cal_f en todas y cada una de sus llamadas, aun cuando implicaba
c        redundancia.
c***********************************************************************                

        subroutine esfuerzos(max_prop,max_jprop,prop,jprop,max_opcion,
     e          opcion,n_jprop,n_prop,max_var_v,max_var_a,max_param,
     e          param,varint_v_cero,varint_a_cero,max_def,max_fterm,
     e          deftot_cero,deftot,esf_cero,max_act,dtiempo,
     s          varint_v,varint_a,esf,fterm,activas,
     s                  convergencia)
c
       IMPLICIT NONE
c******************************************************************
c        DECLARACION DE VARIABLES MDC
c******************************************************************
c
c        prop= matriz de propiedades del elemento
c        param= parametros de la integracion numerica
c        param(1)=alfa
c        param(2)=beta
c        param(3)=gamma
c        param(4)=mﬂximo nînero de iteraciones
c        param(5)=mﬂximo nîmero de pasos de integracin local
c        param(6)=alfa lçmite
c        param(7:6+nt)=ceros numericos para {R},{T},{G}
c        opcion= opciones del calculo (dinamico o estatico,
c                peque±os o grandes desplazamientos, etc)
c        opcion(1)= "static" o "dynamic"
c        opcion(2)= "small" o "large"
c        dtiempo= incremento de tiempo en el paso
c        varint_v_cero= variables internas con funcin de fluencia al principio del paso
c        varint_a_cero= variables internas sin funcin de fluencia al principio del paso
c        varint_v= variables internas con funcin de fluencia al final del paso
c        varint_a= variables internas sin funcin de fluencia al final del paso
c        nv= nîmero de variables internas con funcin de fluencia
c        na= nîmero de variables internas sin funcin de fluencia
c        ny= nîmero de fuerzas termodinﬂmicas
c        n_prop= nîmero de propiedades del elemento
c        nfi= dimensin de la matriz de deformaciones y de esfuerzos
c        max_prop= tama±o de la matriz de propiedades
c        max_var_v= tama±o de la matriz de var.int. con funcin de fluencia
c        max_var_a= tama±o de la matriz de var.int. sin funcin de fluencia 
c        max_opcion= tama±o de la matriz de opciones
c        max_param= tama±o de la matriz de parametros
c        max_def= tama±o de la matriz de deformaciones
c        esf=esfuerzos al final del paso
c        deftot=deformaciones totales al final del paso
c        esf_cero= esfuerzos al principio del paso
c        deftot_cero= deformaciones totales al principio del paso
c        fterm= Fuerzas termodinamicas asociadas a las variables 
c                    internas al final del paso
c        activas= matriz que indica el estado de todas las variables
c                internas
c        max_act= tama±o de la matriz activas
c
c******************************************************************
c
        integer max_prop,max_jprop,max_var_v,max_var_a,max_opcion
c
        integer max_param,max_def,max_act,max_fterm
c
        Real(kind=8):: prop(max_prop),param(max_param),dtiempo
c
        character(len=8):: opcion(max_opcion)
c
        integer n_jprop,n_prop,jprop(max_jprop)
c
        Real(kind=8):: varint_v(max_var_v),varint_a(max_var_a),
     2               varint_v_cero(max_var_v),varint_a_cero(max_var_a)
c
        Real(kind=8):: deftot_cero(max_def),esf_cero(max_def),
     1         deftot(max_def),esf(max_def)
c
        Real(kind=8):: fterm(max_fterm)
c
        logical convergencia,activas(max_act)
c
c******************************************************************
c                DECLARACIONES LOCALES
c******************************************************************
c
c        max_ite= maximo nîmero de iteraciones permitido
c        tol= error admisible de la solucin
c        i_ite= numero de la iteracion
c        hiper= hipermatriz
c        max_hiper= tama~o de la hipermatriz
c        rtg= -({R},{T},{G})= funciones que deben igualarse a cero
c        convergencia= variable logica que indica la convergencia o no
c                del proceso de resolucion del sistema de ecuaciones
c        cambio= variable logica que indica un cambio de activas a pasivas
c                o viceversa en la subrutina verif y no en residual
c
c******************************************************************
c
        integer nfi,nv,na,nt
        integer max_hiper,max_f,max_r,max_g
        parameter(max_hiper=11)
        parameter(max_r=3)
        parameter(max_f=6)
        parameter(max_g=2)
        integer max_ite,i_ite,i
        integer max_ite_cambio,ite_cambio
        logical cambio
        Real(kind=8):: hiper(max_hiper,max_hiper),rtg(max_hiper)
        Real(kind=8):: f(max_f),r(max_r),g(max_g)
        logical auxconv
c
c******************************************************************
c
        auxconv=.true.
        nfi=jprop(4)
        nv=jprop(5)
        na=jprop(6)
        nt=nfi+nv+na
c
c***************************Actualizacion de variables ************
        do i=1,nv
               activas(i)=.false.
        end do
c
        do i=1,nfi
               esf(i)=esf_cero(i)
        end do
c
        do i=1,nv
               varint_v(i)=varint_v_cero(i)
        end do
c
        do i=1,na
               varint_a(i)=varint_a_cero(i)
        end do
c
        cambio=.false.
c
c******************************************************************
c        LLAMADA A LA SUBRUTINA CAL_F QUE CALCULA EL VALOR DE
c        {f}  PARA CADA VARIABLE INTERNA V, DONDE {f} ES LA MATRIZ DE
c        FUNCIONES DE FLUENCIA
c******************************************************************
c        
        call cal_f(max_prop,prop,max_opcion,opcion,
     e             max_var_v,max_var_a,n_prop,
     e             max_param,param,varint_v_cero,varint_a_cero,
     e             max_def,max_jprop,jprop,
     e             deftot_cero,deftot,esf_cero,
     e             varint_v,varint_a,esf,max_f,
     s             f,auxconv)

        if (.not. auxconv) then
        convergencia=.false.
        return
        endif
c
c
c******************************************************************
c        LLAMADA A LA SUBRUTINA CAL_R QUE CALCULA EL  VALOR  DE
c        {R}   DONDE  {R} ES LA LEY DE ESTADO QUE DEFINE LA RELACION
c        ENTRE:
c        esf,deftot,varint_v y varint_a
c******************************************************************
c
        call cal_r(max_prop,prop,max_opcion,opcion,
     e           max_jprop,jprop,n_prop,dtiempo,
     e           max_var_v,max_var_a,max_param,param,varint_v_cero,
     e           varint_a_cero,max_def,
     e           deftot_cero,deftot,esf_cero,
     e           varint_v,varint_a,esf,max_r,
     s           r)
c
c******************************************************************
c        LLAMADA A LA SUBRUTINA CAL_G QUE CALCULA EL  VALOR  DE
c        {G}   DONDE  {G} DEFINE LA LEY DE DAN~O
c******************************************************************

        call cal_g(max_prop,prop,max_opcion,opcion,
     e             max_jprop,jprop,n_prop,
     e             max_var_v,max_var_a,max_param,param,
     e             varint_v,varint_a,
     e             varint_v_cero,varint_a_cero,max_def,
     e             deftot_cero,deftot,esf_cero,
     e             esf,max_g,
     s             g)
c
c******************************************************************
c        LLAMADA A LA SUBRUTINA QUE CALCULA Y ENSAMBLA EL VECTOR RTG
c******************************************************************
c
        call cal_rtg(max_prop,prop,max_opcion,opcion,
     e               max_jprop,jprop,n_prop,dtiempo,
     e               max_var_v,max_param,param,
     e               varint_v_cero,max_def,
     e               deftot_cero,deftot,esf_cero,
     e               varint_v,esf,max_hiper,
     e               max_act,max_f,max_r,max_g,f,r,g,activas,
     s               rtg)
c
c*****************************************************************
c
        ite_cambio=0
c
        max_ite_cambio=param(4)
c
        do while ((cambio.and.ite_cambio.lt.max_ite_cambio).or.
     1                    (ite_cambio.eq.0))
c
           ite_cambio=ite_cambio+1
c
           convergencia=.false.
c
           i_ite=0
c
           max_ite=param(4)
c
           do while (.not.convergencia.and.i_ite.le.max_ite)
c
             i_ite=i_ite+1
c
             call cal_hiper(max_prop,prop,max_opcion,opcion,
     e                       max_jprop,jprop,n_prop,
     e                       max_var_v,max_var_a,max_param,param,
     e                       varint_v_cero,varint_a_cero,
     e                       max_def,deftot_cero,deftot,esf_cero,
     e                       varint_v,varint_a,esf,max_hiper,
     e                       activas,max_act,
     s                       hiper,auxconv)
c
             if (.not. auxconv) then
                convergencia=.false.
                return
             endif
c****************************************************************************
c         RESOLUCION DEL SISTEMA DE ECUACIONES
c****************************************************************************
c
             call rfssedp(hiper,rtg,nt)
c
c*****************************************************************
c        LLAMADA A LA SUBRUTINA QUE ACTUALIZA LOS ESFUERZOS Y LAS
c        VARIABLES INTERNAS
c*****************************************************************
c
             call actualizacion(max_prop,prop,max_jprop,jprop,
     e                 n_prop,max_var_v,max_var_a,max_param,param,
     e                 max_def,max_hiper,rtg,
     s                 esf,varint_v,varint_a)
c
c******************************************************************
c        LLAMADA A LA SUBRUTINA CAL_F QUE CALCULA EL VALOR DE
c        {f}  PARA CADA VARIABLE INTERNA V, DONDE {f} ES LA MATRIZ DE
c        FUNCIONES DE FLUENCIA
c******************************************************************
c        
             call cal_f(max_prop,prop,max_opcion,opcion,
     e                  max_var_v,max_var_a,n_prop,
     e                  max_param,param,varint_v_cero,varint_a_cero,
     e                  max_def,max_jprop,jprop,
     e                  deftot_cero,deftot,esf_cero,
     e                  varint_v,varint_a,esf,max_f,
     s                  f,auxconv)
c
c
c
c******************************************************************
c        LLAMADA A LA SUBRUTINA CAL_R3D QUE CALCULA EL  VALOR  DE
c        {R}   DONDE  {R} ES LA LEY DE ESTADO QUE DEFINE LA RELACION
c        ENTRE:
c        esf,deftot,varint_v,varint_a
c******************************************************************
c
             call cal_r(max_prop,prop,max_opcion,opcion,
     e                max_jprop,jprop,n_prop,dtiempo,
     e                max_var_v,max_var_a,max_param,param,varint_v_cero,
     e                varint_a_cero,max_def,
     e                deftot_cero,deftot,esf_cero,
     e                varint_v,varint_a,esf,max_r,
     s                r)
c
c******************************************************************
c        LLAMADA A LA SUBRUTINA CAL_G QUE CALCULA EL  VALOR  DE
c        {G}   DONDE  {G} DEFINE LA LEY DE DAN~O
c******************************************************************
c
             call cal_g(max_prop,prop,max_opcion,opcion,
     e                  max_jprop,jprop,n_prop,
     e                  max_var_v,max_var_a,max_param,param,
     e                  varint_v,varint_a,
     e                  varint_v_cero,varint_a_cero,max_def,
     e                  deftot_cero,deftot,esf_cero,
     e                  esf,max_g,
     s                  g)
c
c
c******************************************************************
c        LLAMADA A LA SUBRUTINA QUE CALCULA Y ENSAMBLA EL VECTOR RV
c******************************************************************
c
             call cal_rtg(max_prop,prop,max_opcion,opcion,
     e                  max_jprop,jprop,n_prop,dtiempo,
     e                  max_var_v,max_param,param,
     e                  varint_v_cero,max_def,
     e                  deftot_cero,deftot,esf_cero,
     e                  varint_v,esf,max_hiper,
     e                  max_act,max_f,max_r,max_g,f,r,g,activas,
     s                  rtg)
c
c
c*****************************************************************
c        LLAMADA A LA SUBRUTINA QUE VERIFICA SI RTG(I)>0 Y EN
c        CONSECUENCIA LA CONVERGENCIA        
c*****************************************************************
c
             call conver(rtg,param,max_param,max_hiper,max_jprop,
     e                  jprop,max_def,cambio,
     s                  convergencia)
c
             if (.not. auxconv) then
                convergencia=.false.
                i_ite=max_ite+1
                ite_cambio=max_ite_cambio+1
             endif
c
           end do
c
           if (i_ite.ge.max_ite) then
                    return
           end if
c
c******************************************************************
c        LLAMADA A LA SUBRUTINA QUE CHEQUEA CUALES VARIABLES INTERNAS
c        ESTAN ACTIVAS SI (VI - VI_0) < TOL
c******************************************************************
c
           call verif(max_prop,prop,max_opcion,opcion,
     e           max_jprop,jprop,n_prop,dtiempo,
     e           max_var_v,max_var_a,max_param,param,varint_v_cero,
     e           varint_a_cero,max_def,deftot_cero,deftot,esf_cero,
     e           varint_v,varint_a,esf,max_act,max_f,f,
     s           activas,cambio)
c
c******************************************************************
c        LLAMADA A LA SUBRUTINA QUE CHEQUEA CUALES VARIABLES INTERNAS
c        ESTAN ACTIVAS SI H(I) > 0
c******************************************************************
c
           call activ(max_prop,prop,
     e            max_jprop,jprop,n_prop,
     e            max_var_v,max_param,param,max_def,
     e            max_act,max_f,f,
     s            activas,cambio)
c
c
       end do
c
       if (ite_cambio.ge.max_ite_cambio) then
                convergencia=.false.
                return
       end if
c
       call cal_fterm(max_prop,prop,max_opcion,opcion,
     e              max_fterm,max_jprop,jprop,n_prop,dtiempo,
     e              max_var_v,max_var_a,max_param,param,varint_v_cero,
     e              max_def,deftot_cero,deftot,esf_cero,
     e              varint_v,varint_a,esf,
     s              fterm)
c
       return
       end
c
c******************************************************************
c
c       SUBRUTINA  CAL_f
c       CALCULO DE {F}
c
c       DONDE {varint_v}=(fpi,fpj,di+,dj+,di-,dj-)          
c             {f}=(f1,f2,f3,f4,f5,f6)
c  
c       f_pos=  M - (1- d+)*alfa+*c+*fp - (1 - d+)*((1 - alfa+)* c+ *p + My+)
c
c       f_neg= - M + (1-d-)* alfa- * c- *fp) - (1 - d-)*((1 - alfa-)* c- *p + My-)
c
c       f1= funcin de fluencia asociada a la variable fpi 
c       f2= funcin inelﬂstica asociada a la variable fpj      
c
c       Da±os Positivos:
c   
c
c       g+  =  dd+  =  0                                si    G+ < Gcr +                     
c
c       g+ =(deltaG+)*(G+/R+)**(nf+)-((dR+/dd)*dd+) = 0           
c
c     
c       donde:
c
c       dd+= incremento del danio positivo
c
c                  deltaG+     si      deltaG+ > 0
c       deltaG+ = 
c                  0           si      deltaG+ < 0
c
c
c       Dan-os Negativos:
c
c       g-  =  dd-  =  0                                si    G- < Gcr -                     
c
c       g- =(deltaG-)*(G-/R-)**(nf-)-((dR-/dd)*dd-) = 0          
c
c       donde:
c
c       dd-= incremento del danio negativo
c
c                  deltaG-     si      deltaG- > 0
c       deltaG- = 
c                  0           si      deltaG- < 0
c
c         
c
c****************************************************************************
c       Escrita por Julio Flrez-Lpez y Marça Eugenia Marante
c       (mayo 2002)
c***************************************************************************
c
        subroutine cal_f(max_prop,prop,max_opcion,opcion,
     e                 max_var_v,max_var_a,n_prop,
     e                 max_param,param,varint_v_cero,varint_a_cero,
     e                 max_def,max_jprop,jprop,
     e                       deftot_cero,deftot,esf_cero,   
     e                 varint_v,varint_a,esf,max_f,
     s                 f,auxconv)
c       
       IMPLICIT NONE
c******************************************************************
c               DECLARACION DE VARIABLES NUCLEO
c******************************************************************
c
c       prop= matriz de propiedades del elemento
c       param= parametros de la integracion numerica
c       param(1)=alfa
c       param(2)=beta
c       param(3)=gamma
c       param(4)=maximo numero de iteraciones
c       param(5)=maximo numero de pasos de integracion=PNEWDT
c       param(6)=alfa lçmite
c         param(7:6+nt)= ceros numericos
c       opcion= opciones del calculo (dinamico o estatico,
c               peque~os o grandes desplazamientos, etc)
c       opcion(1)= "static" o "dynamic"
c       opcion(2)= "small" o "large"
c       dtiempo= incremento de tiempo en el paso
c       varint_a_cero= variables internas "Ai" al principio del paso
c       varint_v_cero= variables internas "vi" al principio del paso
c         varint_a= variables internas "Ai" al final del paso
c         varint_v= variables internas "vi" al final del paso
c       nv= numero de variables internas con funcion de fluencia
c         na= numero de variables internas sin funcion de fluencia
c       ny= numero de fuerzas termodinamicas
c       n_prop= numero de propiedades del elemento
c       max_prop= tama~o de la matriz de propiedades
c       max_var_v= tama~o de la matriz de variables internas "Vi"
c         max_var_a= tama~o de la matriz de variables internas "Ai"
c       max_opcion= tama~o de la matriz de opciones
c       max_param= tama~o de la matriz de parametros
c       max_def= tama~o de la matriz de deformaciones
c       esf=esfuerzos al final del paso
c       deftot=deformaciones totales al final del paso
c       esf_cero= esfuerzos al principio del paso
c       deftot_cero= deformaciones totales al principio del paso
c       activas= matriz que indica el estado de las variables internas
c
c
c*********************************************************************
c
        integer max_prop,max_jprop,jprop(max_jprop)
        integer max_var_v,max_var_a,max_opcion
        integer max_param,max_f,max_def
        integer n_prop
c
        character(len=8):: opcion(max_opcion)
c
        Real(kind=8):: prop(max_prop),param(max_param)
        Real(kind=8):: deftot_cero(max_def),deftot(max_def)
        Real(kind=8):: varint_a(max_var_a),varint_v(max_var_v)
        Real(kind=8):: varint_a_cero(max_var_a),varint_v_cero(max_var_v)
        Real(kind=8):: esf(max_def),esf_cero(max_def)
        Real(kind=8):: f(max_f)
        logical auxconv
c
c*********************************************************************
c               DECLARACIONES LOCALES
c*********************************************************************
c
c       fpi= rotacion plastica en el nudo i
c       fpj= rotacion plastica en el nudo j
c       di_pos= da~o positivo en el nudo i
c       dj_pos= da~o positivo en el nudo j
c       di_neg= da~o negativo en el nudo i
c       dj_neg= da~o negativo en el nudo j
c       pi= valor absoluto de la deformacion plastica maxima en i
c       pj= valor absoluto de la deformacion plastica maxima en j
c       ei= modulo de elasticidad* momento de inercia=   prop(10)
c       ea= modulo de elasticidad* area=                 prop(11)
c       c_pos= coeficiente del criterio de plasticidad=  prop(12)
c       c_neg= coeficiente del criterio de plasticidad=  prop(13)
c       my_pos= coeficiente del criterio de plasticidad= prop(14)
c       my_neg= coeficiente del criterio de plasticidad= prop(15)
c       gcr_pos=coeficiente del criterio de plasticidad= prop(16)
c       gcr_neg=coeficiente del criterio de plasticidad= prop(17)
c       q_pos=coeficiente del criterio de plasticidad=   prop(18)
c       q_neg=coeficiente del criterio de plasticidad=   prop(19)
c       alfa_pos=coeficiente del criterio de plasticidad=prop(20)
c       alfa_neg=coeficiente del criterio de plasticidad=prop(21)
c       beta_pos=1-alfa_pos
c       beta_neg=1-alfa_neg
c       long= longitud del elemento
c       mi= momento en el extremo i
c       mj= momento en el extremo j
c       n= fuerza axial
c       fi= rotacion total en el extremo i
c       fj= rotacion total en el extremo j
c       fa= alargamiento de la cuerda
c       nf= el valor de N para la fatiga de bajo ciclaje.
c
c*********************************************************************
c
c
        Real(kind=8):: fpi,fpj,pi,pj,mi,mj,n,fi,fj,fa
        Real(kind=8):: f0,ei,ea,long
        Real(kind=8):: c_pos_i,my_pos_i
        Real(kind=8):: c_pos_j,my_pos_j
        Real(kind=8):: c_neg_i,my_neg_i
        Real(kind=8):: c_neg_j,my_neg_j
        Real(kind=8):: gcr_pos_i,q_pos_i
        Real(kind=8):: gcr_pos_j,q_pos_j
        Real(kind=8):: gcr_neg_i,q_neg_i
        Real(kind=8):: gcr_neg_j,q_neg_j
        Real(kind=8):: di_pos,dj_pos,di_neg,dj_neg
        Real(kind=8):: du_pos_i,du_neg_i,du_pos_j,du_neg_j
        Real(kind=8):: alfa_pos_i,alfa_neg_i,beta_pos_i,beta_neg_i
        Real(kind=8):: alfa_pos_j,alfa_neg_j,beta_pos_j,beta_neg_j
        Real(kind=8):: mi_cero,mj_cero,n_cero
        Real(kind=8):: di_pos_cero,dj_pos_cero,di_neg_cero,dj_neg_cero
        Real(kind=8):: fpi_cero,fpj_cero,pi_cero,pj_cero
        Real(kind=8):: cerog,cero1,cero2
        Real(kind=8):: factor
        integer nfi,nv,na,i
c
c*********************************************************************
c               TRADUCCION NUCLEO --> MDC
c*********************************************************************
c
        long=prop(1)
        ei=prop(4)
        ea=prop(5)
c
        fi=deftot(1)
        fj=deftot(2)
        fa=deftot(3)
c
        mi=esf(1)
        mj=esf(2)
        n=esf(3)
c
        mi_cero=esf_cero(1)
        mj_cero=esf_cero(2)
        n_cero=esf_cero(3)
c
        fpi=varint_v(1)
        fpj=varint_v(2)
        di_pos=varint_v(3)
        dj_pos=varint_v(4)
        di_neg=varint_v(5)
        dj_neg=varint_v(6)
c
        fpi_cero=varint_v_cero(1)
        fpj_cero=varint_v_cero(2)
        di_pos_cero=varint_v_cero(3)
        dj_pos_cero=varint_v_cero(4)
        di_neg_cero=varint_v_cero(5)
        dj_neg_cero=varint_v_cero(6)
c
        pi=varint_a(1)
        pj=varint_a(2)
c
        pi_cero=varint_a_cero(1)
        pj_cero=varint_a_cero(2)
c
        c_pos_i=prop(6)
        c_neg_i=prop(14)
        my_pos_i=prop(8)
        my_neg_i=prop(16)
c
        nfi=jprop(4)
        nv=jprop(5)
        na=jprop(6)
c
        factor=0.001
c
c********************** alfa ===> Cinematico *******
        alfa_pos_i=prop(55)
        alfa_neg_i=prop(55)
c*********************  beta ====> Isotropo
        beta_pos_i=1.-alfa_pos_i
        beta_neg_i=1.-alfa_neg_i
c*********************
        c_pos_j=prop(10)
        c_neg_j=prop(18)
        my_pos_j=prop(12)
        my_neg_j=prop(20)
c**********************
        alfa_pos_j=prop(55)
        alfa_neg_j=prop(55)
        beta_pos_j=1.-alfa_pos_j
        beta_neg_j=1.-alfa_neg_j
c
        gcr_pos_i=prop(9)
        gcr_neg_i=prop(17)
        q_pos_i=prop(7)
        q_neg_i=prop(15)
c
        gcr_pos_j=prop(13)
        gcr_neg_j=prop(21)
        q_pos_j=prop(11)
        q_neg_j=prop(19)
c
        du_pos_i=0.63
        du_neg_i=0.63
        du_pos_j=0.63
        du_neg_j=0.63
c
       f0=long/(3.*ei)
c
c*********************************************************************
c               CALCULO DE LAS FUNCIONES DE FLUENCIA
c*********************************************************************
c
c
c      *****************************************************
c      ****************     EN EL NODO i    ****************
c      *****************************************************
c
         call cal_plast(mi,di_pos,di_neg,alfa_pos_i,
     e                  alfa_neg_i,c_pos_i,c_neg_i,fpi,pi,
     e                  beta_pos_i,beta_neg_i,
     e                  my_pos_i,my_neg_i,
     s                  f(1))
c
c
c      *****************************************************
c      ****************     EN EL NODO j    ****************
c      *****************************************************
c
         call cal_plast(mj,dj_pos,dj_neg,alfa_pos_j,
     e                  alfa_neg_j,c_pos_j,c_neg_j,fpj,pj,
     e                  beta_pos_j,beta_neg_j,
     e                  my_pos_j,my_neg_j,
     s                  f(2))
c
c
c*********************************************************************
c               CALCULO DE LAS FUNCIONES DE DA“O
c*********************************************************************
c
c
c      *****************************************************
c      ****   EN EL NODO i, Da±o positivo y negativo    ****
c      *****************************************************

         call cal_dano(f0,mi,di_pos,di_neg,gcr_pos_i,gcr_neg_i,
     e             q_pos_i,q_neg_i,factor,
     e                     du_pos_i,du_neg_i,
     e             mi_cero,di_pos_cero,di_neg_cero,
     s             f(3),f(5),cero1,auxconv)


c         call cal_dano(f0,mi,di_pos,di_neg,gcr_pos_i,gcr_neg_i,
c     e             q_pos_i,q_neg_i,factor,
c     e                     du_pos_i,du_neg_i,
c     e             mi_cero,di_pos_cero,di_neg_cero,
c     s             f(3),f(5),cero1)

c
c
c
c      *****************************************************
c      ****   EN EL NODO j, Da±o positivo y negativo    ****
c      *****************************************************

         call cal_dano(f0,mj,dj_pos,dj_neg,gcr_pos_j,gcr_neg_j,
     e             q_pos_j,q_neg_j,factor,
     e                     du_pos_j,du_neg_j,
     e             mj_cero,dj_pos_cero,dj_neg_cero,
     s             f(4),f(6),cero2,auxconv)


c         call cal_dano(f0,mj,dj_pos,dj_neg,gcr_pos_j,gcr_neg_j,
c     e             q_pos_j,q_neg_j,factor,
c     e                     du_pos_j,du_neg_j,
c     e             mj_cero,dj_pos_cero,dj_neg_cero,
c     s             f(4),f(6),cero2)

c
c*********************************************************************
c         CALCULO DEL CERO NUMERICO PARA LAS FUNCIONES DE DANO
c*********************************************************************
c
c
         if(cero1.gt.cero2) then 
            cerog=cero1
         else
            cerog=cero2
         endif
c
         do i=nfi+nv-3,nfi+nv
            param(i+6)=cerog*factor
         end do
c
        return
        end
c
c*********************************************************************
c           SUBRUTINA CAL_PLAST
c           CALCULO DE LAS FUNCIONES DE FLUENCIA EN 
c           ELEMENTOS DE CONCRETO ARMADO.
c
c                  Modificada Nayive Jaramillo (Enero 2002) 
c*********************************************************************
c
c
         subroutine cal_plast(m,d_pos,d_neg,alfa_pos,
     e                 alfa_neg,c_pos,c_neg,fp,p,
     e                 beta_pos,beta_neg,my_pos,my_neg,
     s                 ff)
c
c
       IMPLICIT NONE
c*********************************************************************
c               DECLARACIONES
c*********************************************************************
c
c       m= momento en el extremo
c       d_pos= da±o positivo en el nudo
c       d_neg= da±o negativo en el nudo
c       alfa_pos=coeficiente del criterio de plasticidad=prop(20)
c       alfa_neg=coeficiente del criterio de plasticidad=prop(21)
c       c_pos= coeficiente del criterio de plasticidad=  prop(12)
c       c_neg= coeficiente del criterio de plasticidad=  prop(13)
c       fp= rotacion plastica en el nudo
c       p= valor absoluto de la deformacion plastica maxima
c       beta_pos=1-alfa_pos
c       beta_neg=1-alfa_neg
c       my_pos= coeficiente del criterio de plasticidad= prop(14)
c       my_neg= coeficiente del criterio de plasticidad= prop(15)
c       ff_pos= funcion de fluencia positiva en el nodo
c       ff_neg= funcion de fluencia negativo en el nodo
c
c*********************************************************************
c
        Real(kind=8):: fp,m,p
        Real(kind=8):: c_pos,my_pos,c_neg,my_neg
        Real(kind=8):: d_pos,d_neg
        Real(kind=8):: alfa_pos,alfa_neg,beta_pos,beta_neg
        Real(kind=8):: ff,valor_pos,valor_neg
c
c********************************************************************
c
c
        valor_pos=m/(1.-d_pos)/(1.-d_neg)-alfa_pos*c_pos*fp
c
        valor_neg=-m/(1.-d_pos)/(1.-d_neg)+alfa_neg*c_neg*fp
c
c
        if(valor_pos.ge.valor_neg)then
c
                ff=m-(1.-d_pos)*alfa_pos*c_pos*fp-
     1             (1.-d_pos)*
     2          (beta_pos*c_pos*p+my_pos)
c
        else
c
                ff=-m+(1.-d_neg)*alfa_neg*c_neg*fp-
     1             (1.-d_neg)*
     2          (beta_neg*c_neg*p+my_neg)
c
        end if
c
        return
c
        end
c
c
c*********************************************************************
c
c           SUBRUTINA CAL_DANO
c           CALCULO DE LAS FUNCIONES DE DA“O, TOMANDO EN CUENTA 
c           LA FATIGA DE BAJO CICLAJE
c
c**********************************************************************
c
         subroutine cal_dano(f0,m,d_pos,d_neg,gcr_pos,gcr_neg,
     e             q_pos,q_neg,factor,
     e             du_pos,du_neg,
     e             m_cero,d_pos_cero,d_neg_cero,
     s             gda1,gda2,cero,auxconv)
c
c
       IMPLICIT NONE
c*********************************************************************
c                 DECLARACIONES
c*********************************************************************
c
c        m= momento en el nodo
c        d_pos= da±o positivo en el nodo
c        d_neg= da±o negativo en el nodo
c        gcr_pos= coeficiente positivo en el nodo del modelo
c        gcr_neg= coeficiente negativo en el nodo del modelo
c        q_pos= coeficiente positivo en el nodo del modelo
c        q_neg= coeficiente negativo en el nodo del modelo
c        du_pos= coeficiente positivo en el nodo del modelo
c        du_neg= coeficiente negativo en el nodo del modelo
c
c        m_cero= momento en el nodo del paso anterior
c        d_pos_cero= da±o positivo en el nodo del paso anterior
c        d_neg_cero= da±o negativo en el nodo del paso anterior
c        gda1= valor de la funcin de da±o positiva en el nodo
c        gda2= valor de la funcin de da±o negativa en el nodo
c        c1,c2,c3,c4=variables auxiliares para realizar el calculo de 
c                    los valores de cero1 y cero2
c
c*********************************************************************
c
c
        Real(kind=8):: m,m_cero,f0,d_pos,d_neg
        Real(kind=8):: d_pos_cero,d_neg_cero,du_pos,du_neg
        Real(kind=8):: gda1,gda2,gcr_pos,gcr_neg,q_pos,q_neg,cero
        Real(kind=8):: factor
c
c
c********************************************************************
c
        Real(kind=8):: no,pos_m,neg_m
        Real(kind=8):: pos_m_cero,neg_m_cero
        Real(kind=8):: kk1_pos,kk2_pos,nf_pos
        Real(kind=8):: kk1_neg,kk2_neg,nf_neg
        Real(kind=8):: G_pos,G_neg,G_pos_cero,G_neg_cero
        Real(kind=8):: delta_G_pos,delta_G_neg
        Real(kind=8):: pos_delta_G_pos,pos_delta_G_neg
        Real(kind=8):: rdpos,derivrpos,rdneg,derivrneg
        Real(kind=8):: c1,c2,c3,c4
        logical auxconv
c
c********************************************************************
c
        no=30.
        auxconv=.true.
c
c**************************************************************
c               CON FATIGA DE BAJO CICLAJE
c**************************************************************
c
        if(m.ge.0.)then
                        pos_m=m
                        neg_m=0.
        else
                        pos_m=0.
                        neg_m=m
        end if
c
        if(m_cero.ge.0.)then
                        pos_m_cero=m_cero
                        neg_m_cero=0.
        else
                        pos_m_cero=0.
                        neg_m_cero=m_cero
        end if
c
c*********************************************************************
c       CALCULO DE "N" PARA LA FATIGA DE BAJO CICLAJE EN EL NODO
c*********************************************************************
c
c       Parte Positiva:

              kk1_pos=(5.+no*(du_pos-1.))/(du_pos*
     1                (du_pos-1.))
              kk2_pos=-kk1_pos-no
              nf_pos=kk1_pos*d_pos*d_pos+
     1               kk2_pos*d_pos+no
c
c
c       Parte Negativa:
c
              kk1_neg=(5.+no*(du_neg-1.))/(du_neg*
     1                (du_neg-1.))
              kk2_neg=-kk1_neg-no
              nf_neg=kk1_neg*d_neg*d_neg+
     1               kk2_neg*d_neg+no

c        if((nf_neg.lt.0.).or.(nf_neg.gt.30.).or.
c     1            (nf_pos.lt.0.).or.(nf_pos.gt.30.)) then
c                auxconv=.false.
c                return
c        endif
c
c.
*********************************************************************
c       CALCULO DE G+  Y  G-
c*********************************************************************
c
              G_pos=pos_m*pos_m*f0/(2.*(1.-d_pos)*(1.-d_pos))
              G_neg=neg_m*neg_m*f0/(2.*(1.-d_neg)*(1.-d_neg))
c
              G_pos_cero=pos_m_cero*pos_m_cero*f0/
     1                  (2.*(1.-d_pos_cero)*(1.-d_pos_cero))
c
              G_neg_cero=neg_m_cero*neg_m_cero*f0/
     1                  (2.*(1.-d_neg_cero)*(1.-d_neg_cero))
c
c*********************************************************************
c
              delta_G_pos=G_pos-G_pos_cero
              delta_G_neg=G_neg-G_neg_cero
c
c
            if(delta_G_pos.gt.0.)then
                 pos_delta_G_pos=delta_G_pos
            else
                 pos_delta_G_pos=0.
            end if
c
            if(delta_G_neg.gt.0.)then
                 pos_delta_G_neg=delta_G_neg
            else
                 pos_delta_G_neg=0.
            end if
c
c                ***********************************************************
c                  Calculo de c1 y de la tasa de disipacion de energia pos
c                ***********************************************************
            if(G_pos.lt.gcr_pos)then
c
                 gda1=d_pos-d_pos_cero
                 c1=factor
            else
           if ((d_pos.gt.1.).or.(d_pos.lt.0)) then
                auxconv=.false.
                return
           endif
                 rdpos=gcr_pos+q_pos*dlog(1.-d_pos)/
     1                 (1.-d_pos)
c
                 derivrpos=-q_pos*(1.-dlog(1.-d_pos))/
     1                          ((1.-d_pos)*(1.-d_pos))

           if ((G_pos/rdpos) .lt. 0.) then
                auxconv=.false.
                return
           endif
 
                gda1=pos_delta_G_pos*(G_pos/rdpos)**(nf_pos)-
     1                 derivrpos*(d_pos-d_pos_cero)          
c
                 c3=pos_delta_G_pos*(G_pos/rdpos)**(nf_pos)
                 c4=derivrpos*(d_pos-d_pos_cero)
c
                 if(dabs(c3).gt.dabs(c4))then
                    c1=dabs(c3)
                 else
                    c1=dabs(c4)
                 end if
c
            end if
c
c                ************************************************************ 
c                  Calculo de c2 y de la tasa de disipacion de energia neg      
c                ************************************************************
c
            if(G_neg.lt.gcr_neg)then
c
                 gda2=d_neg-d_neg_cero
                 c2=factor
            else
           if ((d_neg.gt.1.).or.(d_neg.lt.0.)) then
                auxconv=.false.
                return
           endif
                 rdneg=gcr_neg+q_neg*dlog(1.-d_neg)/(1.-d_neg)
 
                 derivrneg=-q_neg*(1.-dlog(1.-d_neg))/
     1                          ((1.-d_neg)*(1.-d_neg))

            if ((G_neg/rdneg) .lt. 0.) then
                auxconv=.false.
                return
            endif
 
                gda2=pos_delta_G_neg*(G_neg/rdneg)**(nf_neg)-
     1                derivrneg*(d_neg-d_neg_cero)    
c
                 c3=pos_delta_G_neg*(G_neg/rdneg)**(nf_neg)
                 c4=derivrneg*(d_neg-d_neg_cero)
c
                 if(dabs(c3).gt.dabs(c4))then
                    c2=dabs(c3)
                 else
                    c2=dabs(c4)
                 end if
c
            end if
c                *********************************************** 
c                 CALCULO DE cero
c                ***********************************************
c
            if(c1.gt.c2) then
               cero=c1
            else
               cero=c2
            end if
        return
c
        end
c
c*********************************************************************
c
c           SUBRUTINA CAL_G
c           CALCULO DE LA ROTACI‘N PL¬STICA MAXIMA
c
c**********************************************************************
c
c       Escrita por Marça Eugenia Marante y Julio Flrez-Lpez
c       (Mayo 2002)
c
c******************************************************************
c
        subroutine cal_g(max_prop,prop,max_opcion,opcion,
     e                        max_jprop,jprop,n_prop,
     e                        max_var_v,max_var_a,max_param,param,
     e                        varint_v,varint_a,
     e                        varint_v_cero,varint_a_cero,max_def,
     e                        deftot_cero,deftot,esf_cero,
     e                        esf,max_g,
     s                           g)            
c******************************************************************
c               DECLARACION DE VARIABLES NUCLEO
c******************************************************************
c
       IMPLICIT NONE
c
        integer max_prop,max_jprop,max_var_v,max_var_a,max_opcion
        integer max_param,max_def,max_g
        integer n_prop,jprop(max_jprop)
c
        character(len=8):: opcion(max_opcion)
c
        Real(kind=8):: prop(max_prop),param(max_param)
        Real(kind=8):: varint_v(max_var_v),varint_v_cero(max_var_v)
        Real(kind=8):: varint_a(max_var_a),varint_a_cero(max_var_a)
        Real(kind=8):: deftot_cero(max_def),esf_cero(max_def)
             Real(kind=8):: deftot(max_def),esf(max_def)
        Real(kind=8):: g(max_g)
c
c
c*********************************************************************
c               DECLARACIONES LOCALES
c*********************************************************************
c
c       fpi= rotacion plastica en el nudo i
c       fpj= rotacion plastica en el nudo j
c       pi= valor absoluto de la deformacion plastica maxima en i
c       pj= valor absoluto de la deformacion plastica maxima en j
c
c*********************************************************************
c
        Real(kind=8):: fpi,fpj,pi,pj
          Real(kind=8):: pi_cero,pj_cero
c
c*********************************************************************
c               TRADUCCION NUCLEO --> MDC
c*********************************************************************
c
        fpi=varint_v(1)
        fpj=varint_v(2)
c
          pi=varint_a(1)
        pj=varint_a(2)
c
          pi_cero=varint_a_cero(1)
          pj_cero=varint_a_cero(2)
c
c        
          if(dabs(fpi).ge.pi_cero)then
                        g(1)=pi-dabs(fpi)
          else
                        g(1)=pi-pi_cero
          end if
c
          if(dabs(fpj).ge.pj_cero)then
                        g(2)=pj-dabs(fpj)
          else
                        g(2)=pj-pj_cero
          end if
c
        return
        end
c

c******************************************************************
c
c       SUBRUTINA  CAL_R
c       CALCULO DE {R}
c
c       {R}={f-fp}-[F(D+)]{M+}-[F(D-)]{M-}
c
c       J. FLOREZ LOPEZ
c       A. RAMIREZ
c       (Nov./94)
c
c
c
c******************************************************************
c
        subroutine cal_r(max_prop,prop,max_opcion,opcion,
     e             max_jprop,jprop,n_prop,dtiempo,
     e             max_var_v,max_var_a,max_param,param,varint_v_cero,
     e             varint_a_cero,max_def,
     e             deftot_cero,deftot,esf_cero,
     e             varint_v,varint_a,esf,max_r,
     s               r)
c
       IMPLICIT NONE
c******************************************************************
c               DECLARACION DE VARIABLES NUCLEO
c******************************************************************
c
c       prop= matriz de propiedades del elemento
c       param= parametros de la integracion numerica
c       param(1)=alfa
c       param(2)=beta
c       param(3)=gamma
c       param(4)=maximo numero de iteraciones
c       param(5)=maximo numero de pasos de integ.
c       param(6)=alfa limite
c
c       opcion= opciones del calculo (dinamico o estatico,
c               peque~os o grandes desplazamientos, etc)
c       opcion(1)= "static" o "dynamic"
c       opcion(2)= "small" o "large"
c       dtiempo= incremento de tiempo en el paso
c       varint_v_cero= variables internas "vi" al principio del paso
c       varint_v= variables internas "vi" al final del paso
c       n_prop= numero de propiedades del elemento
c          varint_a_cero= variables internas "ai" al principio del paso
c       varint_a= variables internas "ai" al final del paso
c       max_prop= tama~o de la matriz de propiedades
c       max_var_a= tama~o de la matriz de variables internas "ai"
c       max_var_v= tama~o de la matriz de variables internas "vi"
c       max_opcion= tama~o de la matriz de opciones
c       max_param= tama~o de la matriz de parametros
c       max_def= tama~o de la matriz de deformaciones
c       esf=esfuerzos al final del paso
c       deftot=deformaciones totales al final del paso
c       esf_cero= esfuerzos al principio del paso
c       deftot_cero= deformaciones totales al principio del paso
c       rgt=({R},{V})
c       activas= matriz que indica el estado de las variables internas
c
c*********************************************************************
c
        integer max_prop,max_opcion,max_param
        integer max_def,max_r
        integer n_prop,max_jprop,jprop 
          integer max_var_v,max_var_a
c
        Real(kind=8):: varint_v_cero(max_var_v),varint_a_cero(max_var_a)
        character(len=8):: opcion(max_opcion)
c
        Real(kind=8):: varint_a(max_var_a),varint_v(max_var_v)
        Real(kind=8):: deftot_cero(max_def),esf_cero(max_def)
        Real(kind=8):: deftot(max_def),esf(max_def)
        Real(kind=8):: r(max_r)
        Real(kind=8):: prop(max_prop),param(max_param),dtiempo
c
c*********************************************************************
c               DECLARACIONES LOCALES
c*********************************************************************
c
c       fpi= rotacion plastica en el nudo i
c       fpj= rotacion plastica en el nudo j
c       di_pos= da~o positivo en el nudo i
c       dj_pos= da~o positivo en el nudo j
c       di_neg= da~o negativo en el nudo i
c       dj_neg= da~o negativo en el nudo j
c       ei= modulo de elasticidad* momento de inercia=prop(10)
c       ea= modulo de elasticidad* area=prop(11)
c       long= longitud del elemento
c       mi= momento en el extremo i
c       mj= momento en el extremo j
c       n= fuerza axial
c       fi= rotacion total en el extremo i
c       fj= rotacion total en el extremo j
c       fa= alargamiento de la cuerda
c       f11_pos,f12_pos,f22_pos,f33= elementos de la matriz de 
c                                    flexibilidad con da~o positivo
c       f11_neg,f12_neg,f22_neg,f33= elementos de la matriz de 
c                                    flexibilidad con da~o negativo
c
c*********************************************************************
c
        Real(kind=8):: fpi,fpj,ei,ea,long,mi,mj,n,fi,fj,fa
        Real(kind=8):: f11_pos,f12_pos,f22_pos
        Real(kind=8):: f11_neg,f12_neg,f22_neg,f33
        Real(kind=8):: di_pos,dj_pos,di_neg,dj_neg
        Real(kind=8):: pos_mi,neg_mi,pos_mj,neg_mj
c
c*********************************************************************
c               TRADUCCION NUCLEO --> MOD1
c*********************************************************************
c
        fi=deftot(1)
        fj=deftot(2)
        fa=deftot(3)
c
        mi=esf(1)
        mj=esf(2)
        n=esf(3)
c
        fpi=varint_v(1)
        fpj=varint_v(2)
        di_pos=varint_v(3)
        dj_pos=varint_v(4)
        di_neg=varint_v(5)
        dj_neg=varint_v(6)
c
        long=prop(1)
        ei=prop(4)
        ea=prop(5)
c
c*********************************************************************
c               CALCULO DE LAS MATRICES DE FLEXIBILIDAD
c*********************************************************************
c
        f11_pos=long/(3.*ei*(1.-di_pos))
        f12_pos=-long/(6.*ei)
        f22_pos=long/(3.*ei*(1.-dj_pos))
c
        f11_neg=long/(3.*ei*(1.-di_neg))
        f12_neg=-long/(6.*ei)
        f22_neg=long/(3.*ei*(1.-dj_neg))
c
        f33=long/ea
c
c*********************************************************************
c               CALCULO DE LAS PARTES POSITIVA Y NEGATIVA DE LOS
c               MOMENTOS
c*********************************************************************
c
        if(mi.ge.0)then
                pos_mi=mi
                neg_mi=0.
        else
                pos_mi=0.
                neg_mi=mi
        end if
c
        if(mj.ge.0.)then
                pos_mj=mj
                neg_mj=0.
        else
                pos_mj=0.
                neg_mj=mj
        end if
c
c*********************************************************************
c               CALCULO DE {R}
c*********************************************************************
c
        r(1)=fi-fpi-(f11_pos*pos_mi+f12_pos*pos_mj)
     1                  -(f11_neg*neg_mi+f12_neg*neg_mj)
c
        r(2)=fj-fpj-(f12_pos*pos_mi+f22_pos*pos_mj)
     1                  -(f12_neg*neg_mi+f22_neg*neg_mj)
c
        r(3)=fa-f33*n
c
        return
        end
c
c******************************************************************
c
c        SUBRUTINA  CAL_RTG
c
c        ESCRITA POR Marça Eugenia Marante y Julio Flrez Lpez
c
c        (Marzo 2001)
c     Modificada Nayive Jaramillo (Enero 2002)
c******************************************************************
c
        subroutine cal_rtg(max_prop,prop,max_opcion,opcion,
     e                   max_jprop,jprop,n_prop,dtiempo,
     e                   max_var_v,max_param,param,
     e                   varint_v_cero,max_def,
     e                   deftot_cero,deftot,esf_cero,
     e                   varint_v,esf,max_hiper,
     e                   max_act,max_f,max_r,max_g,f,r,g,activas,
     s                    rtg)
c
       IMPLICIT NONE
c******************************************************************
c        DECLARACION DE VARIABLES MDC
c******************************************************************
c
c        prop= matriz de propiedades del elemento
c        param= parametros de la integracion numerica
c        param(1)=alfa
c        param(2)=beta
c        param(3)=gamma
c        param(4)=maximo numero de iteraciones
c        param(5)=maximo numero de pasos de integracin local
c        param(6:nfi+nv+na)=ceros numâricos de {R},{T},{G}
c        opcion= opciones del calculo (dinamico o estatico,
c                peque~os o grandes desplazamientos, etc)
c        opcion(1)= "static" o "dynamic"
c        opcion(2)= "small" o "large"
c        dtiempo= incremento de tiempo en el paso
c        varint_v_cero= variables internas con funcin de fluencia al principio del paso
c        varint_v = variables internas con funcin de fluencia al final del paso
c        varint_a_cero= variables internas sin funcin de fluencia al principio del paso
c        varint_a = variables internas sin funcin de fluencia al final del paso
c        nv= nîmero de variables internas con funcin de fluencia
c        na= nîmero de variables internas sin funcin de fluencia
c        nfi= dimensin de la matriz de deformaciones y de esfuerzos
c        n_prop= numero de propiedades del elemento
c        max_prop= tama~o de la matriz de propiedades
c        max_var_v= tama~o de la matriz de variables internas con funcin de fluencia
c        max_opcion= tama~o de la matriz de opciones
c        max_param= tama~o de la matriz de parametros
c        max_def= tama~o de la matriz de deformaciones
c        esf=esfuerzos al final del paso
c        deftot=deformaciones totales al final del paso
c        esf_cero= esfuerzos al principio del paso
c        deftot_cero= deformaciones totales al principio del paso
c
c******************************************************************
c
        integer max_prop,max_jprop,max_var_v,max_opcion,max_param
c
        integer max_def,max_hiper,max_act,max_f,max_r,max_g
c
        Real(kind=8):: prop(max_prop),param(max_param),
     2                dtiempo,f(max_f),r(max_r),g(max_g)
c
        character(len=8):: opcion(max_opcion)
c
        integer n_prop,jprop(max_jprop)
c
        Real(kind=8):: varint_v(max_var_v),
     2                        varint_v_cero(max_var_v)
c
        Real(kind=8):: deftot_cero(max_def),esf_cero(max_def),
     1                deftot(max_def),esf(max_def)
c
        Real(kind=8)::  rtg(max_hiper)
c
        logical activas(max_act)
c
c******************************************************************
c                DECLARACIONES LOCALES
c******************************************************************
c
        integer max_fi,i,nfi,nv,na
        parameter(max_fi=6)
        Real(kind=8):: t(max_fi)
c
c******************************************************************
c        CALCULO DE Vi
c******************************************************************
c
      nfi=jprop(4)
      nv=jprop(5)
      na=jprop(6)
c
c******************************************************************
c
        do i=1,nv
                if(activas(i))then
                        t(i)=f(i)
                else
                        t(i)=varint_v(i)-varint_v_cero(i)
                end if
      end do
c
c******************************************************************
c        ENSAMBLAJE DE ({R},{T},{G})
c******************************************************************
c
         do i=1,nfi
                rtg(i)=r(i)
         end do
c
         do i=nfi+1,nfi+nv
                rtg(i)=t(i-nfi)
         end do
c        
         do i=nfi+nv+1,nfi+nv+na
                rtg(i)=g(i-nfi-nv)
         end do
c
         return
         end
c
c******************************************************************
c         SUBRUTINA  CAL_HIPER
c        CALCULO DE LA HIPERMATRIZ COMPUESTA POR:
c
c        |@R/@ESF    @R/@INT_V    @R/@INT_A|
c        |                                  |
c        |@T/@ESF    @T/@INT_V    @T/@INT_A|
c       |                                  | 
c        |@G/@ESF    @G/@INT_V    @G/@INT_A|
c        ....
c        DONDE
c
c        {R} representa a la ley de estado
c
c        Ti=  f(i)                                        si la variable interna i esta activa
c        Ti = varint_v(i) - varint_v_cero(i)                si la variable interna i no esta activa 
c
c        f(i) es la funcion de carga de la variable interna Vi
c
c        {G} define la Ley de Da~no
c
c        El simbolo @ representa la derivada parcial, ESF los esfuerzos 
c        
c        INT_V las variables internas con funcin de fluencia
c        INT_A las variables internas sin funcin de fluencia
c
c
c        
c******************************************************************
c
        subroutine cal_hiper(max_prop,prop,max_opcion,opcion,
     e                        max_jprop,jprop,n_prop,
     e                        max_var_v,max_var_a,max_param,param,
     e                        varint_v_cero,varint_a_cero,
     e            max_def,deftot_cero,deftot,esf_cero,
     e                        varint_v,varint_a,esf,max_hiper,
     e                        activas,max_act,
     s                        hiper,auxconv)
c
       IMPLICIT NONE
c******************************************************************
c                DECLARACION DE VARIABLES ESFUERZOS
c******************************************************************
c
c        prop= matriz de propiedades del elemento
c        param= parametros de la integracion numerica
c        param(1)=alfa
c        param(2)=beta
c        param(3)=gamma
c        param(4)=maximo numero de iteraciones
c        param(5)=maximo numero de pasos de integracin local
c        param(7:6+nfi+nv+na)=ceros numâricos de {R},{T},{G}
c        opcion= opciones del calculo (dinamico o estatico,
c                peque~os o grandes desplazamientos, etc)
c        opcion(1)= "static" o "dynamic"
c        opcion(2)= "small" o "large"
c        dtiempo= incremento de tiempo en el paso
c        varint_v_cero= variables internas con funcin de fluencia al principio del paso
c        varint_v = variables internas con funcin de fluencia al final del paso
c        varint_a_cero= variables internas sin funcin de fluencia al principio del paso
c        varint_a = variables internas sin funcin de fluencia al final del paso
c        nv= nîmero de variables internas con funcin de fluencia
c        na= nîmero de variables internas sin funcin de fluencia
c        nfi= dimensin de la matriz de deformaciones y de esfuerzos
c        n_prop= numero de propiedades del elemento
c        max_prop= tama~o de la matriz de propiedades
c        max_var_v= tama~o de la matriz de variables internas con funcin de fluencia
c        max_opcion= tama~o de la matriz de opciones
c        max_param= tama~o de la matriz de parametros
c        max_def= tama~o de la matriz de deformaciones
c        esf=esfuerzos al final del paso
c        deftot=deformaciones totales al final del paso
c        esfuerzo_cero= esfuerzos al principio del paso
c        deftot_cero= deformaciones totales al principio del paso
c        hiper= hipermatriz para la resolucion del sistema de ecuaciones
c               locales por el metodo de Newton
c        max_hiper= tama~o maximo de hiper
c        activas= matriz que indica el estado de las variables internas
c
c******************************************************************
c
        integer max_prop,max_var_v,max_var_a,max_opcion,max_param
        integer max_def,max_hiper,max_act,max_jprop
c
        Real(kind=8):: prop(max_prop),param(max_param)
c
        character(len=8):: opcion(max_opcion)
c
        integer n_prop,jprop(max_jprop)
c
        Real(kind=8):: varint_v(max_var_v),varint_a(max_var_a),
     2           varint_v_cero(max_var_v),varint_a_cero(max_var_a)
c
        Real(kind=8):: deftot_cero(max_def),esf_cero(max_def),
     1                deftot(max_def),esf(max_def)
c
c
        Real(kind=8):: hiper(max_hiper,max_hiper)
c
        logical activas(max_act),auxconv
c
c
c******************************************************************
c                DECLARACIONES LOCALES
c******************************************************************
c
c        deriv_r_esf= matriz de derivadas de {R} con respecto a los
c                        esfuerzos
c
c        deriv_r_intv= matriz de derivadas de {R} con respecto a las
c                        variables internas con funcin de fluencia
c
c        deriv_r_inta= matriz de derivadas de {R} con respecto a las
c                        variables internas sin funcin de fluencia
c
c        deriv_f_esf= matriz de derivadas de {F} con respecto a los
c                        esfuerzos
c
c        deriv_f_intv= matriz de derivadas de {F} con respecto a las
c                        variables internas con funcin de fluencia
c
c        deriv_f_inta= matriz de derivadas de {F} con respecto a las
c                        variables internas sin funcin de fluencia
c
c        deriv_g_esf= matriz de derivadas de {G} con respecto a los
c                        esfuerzos
c
c        deriv_g_intv= matriz de derivadas de {G} con respecto a las
c                        variables internas con funcin de fluencia
c
c        deriv_g_inta= matriz de derivadas de {G} con respecto a las
c                        variables internas sin funcin de fluencia
c
c        max_esf= dimension maxima de deriv_r_esf=max_def
c        max_int_v= dimension maxima de deriv_f_intv=max_var_v
c        max_int_a= dimension maxima de deriv_f_inta=max_var_a
c
c******************************************************************
c
        integer nfi,nv,na,i,j
        integer max_esf,max_int_v,max_int_a
c
        parameter(max_esf=3)
        parameter(max_int_v=6)
        parameter(max_int_a=2)
c
        Real(kind=8):: deriv_r_esf(max_esf,max_esf)
        Real(kind=8):: deriv_r_intv(max_esf,max_int_v)
        Real(kind=8):: deriv_r_inta(max_esf,max_int_a)
        Real(kind=8):: deriv_f_esf(max_int_v,max_esf)
        Real(kind=8):: deriv_f_intv(max_int_v,max_int_v)
        Real(kind=8):: deriv_f_inta(max_int_v,max_int_a)
        Real(kind=8):: deriv_g_esf(max_int_a,max_esf)
        Real(kind=8):: deriv_g_intv(max_int_a,max_int_v)
        Real(kind=8):: deriv_g_inta(max_int_a,max_int_a)
c
c******************************************************************
c        LLAMADA A  LA SUBRUTINA  QUE  CALCULA
c        LAS DERIVADAS DE R CON RESPECTO A ESF
c******************************************************************
c
        call cal_deriv_r_esf(max_prop,prop,max_opcion,opcion,
     e           max_jprop,jprop,n_prop,max_esf,
     e           max_var_v,max_var_a,max_param,param,varint_v_cero,
     e           varint_a_cero,max_def,
     e           deftot_cero,deftot,esf_cero,
     e           varint_v,varint_a,esf,max_hiper,
     s           deriv_r_esf)
c
c******************************************************************
c        LLAMADA A LA SUBRUTINA  QUE CALCULA LAS
c        DERIVADAS  DE  R  CON  RESPECTO A  VI CON FUNCION DE FLUENCIA
c******************************************************************
c
        call cal_deriv_r_intv(max_prop,prop,max_opcion,opcion,
     e         max_jprop,jprop,n_prop,max_esf,max_int_v,
     e         max_var_v,max_var_a,max_param,param,varint_v_cero,
     e         varint_a_cero,max_def,
     e         deftot_cero,deftot,esf_cero,
     e         varint_v,varint_a,esf,max_hiper,
     s         deriv_r_intv)
c
c******************************************************************
c        LLAMADA A LA SUBRUTINA  QUE CALCULA LAS
c        DERIVADAS  DE  R  CON  RESPECTO A  VI SIN FUNCION DE FLUENCIA
c******************************************************************
c
        call cal_deriv_r_inta(max_prop,prop,max_opcion,opcion,
     e         max_jprop,jprop,n_prop,max_esf,max_int_a,
     e         max_var_v,max_var_a,max_param,param,varint_v_cero,
     e         varint_a_cero,max_def,
     e         deftot_cero,deftot,esf_cero,
     e         varint_v,varint_a,esf,max_hiper,
     s         deriv_r_inta)
c
c******************************************************************
c        LLAMADA A LA SUBRUTINA  QUE CALCULA LAS
c        DERIVADAS  DE  f  CON  RESPECTO  A  ESF
c******************************************************************
c        
        call cal_deriv_f_esf(max_prop,prop,max_opcion,opcion,
     e         max_jprop,jprop,n_prop,max_esf,max_int_v,
     e         max_var_v,max_var_a,max_param,param,varint_v_cero,
     e         varint_a_cero,max_def,
     e         deftot_cero,deftot,esf_cero,
     e         varint_v,varint_a,esf,max_hiper,
     s         deriv_f_esf,auxconv)

        if (.not. auxconv) then
                return
        endif
c
c******************************************************************
c        LLAMADA A LA SUBRUTINA  QUE CALCULA LA
c        DERIVADA DE f CON RESPECTO A LAS VI CON FUENCION DE FLUENCIA
c******************************************************************
c        
        call cal_deriv_f_intv(max_prop,prop,max_opcion,opcion,
     e        max_jprop,jprop,n_prop,max_int_v,max_int_a,max_var_v,
     e        max_var_a,max_param,param,varint_v_cero,varint_a_cero,
     e        max_def,deftot_cero,deftot,esf_cero,
     e        varint_v,varint_a,esf,max_hiper,
     s        deriv_f_intv,auxconv)

        if (.not. auxconv) then
                return
        endif

c
c******************************************************************
c        LLAMADA A LA SUBRUTINA  QUE CALCULA LA
c        DERIVADA DE f CON RESPECTO A LAS VI SIN FUENCION DE FLUENCIA
c******************************************************************
c        
        call cal_deriv_f_inta(max_prop,prop,max_opcion,opcion,
     e        max_jprop,jprop,n_prop,max_int_v,max_int_a,max_var_v,
     e        max_var_a,max_param,param,varint_v_cero,varint_a_cero,
     e        max_def,deftot_cero,deftot,esf_cero,
     e        varint_v,varint_a,esf,max_hiper,
     s        deriv_f_inta,auxconv)

        if (.not. auxconv) then
                return
        endif

c
c******************************************************************
c        LLAMADA A LA SUBRUTINA  QUE CALCULA LAS
c        DERIVADAS  DE  g  CON  RESPECTO  A  ESF
c******************************************************************
c        
        call cal_deriv_g_esf(max_prop,prop,max_opcion,opcion,
     e                 max_jprop,jprop,n_prop,max_int_a,max_esf,
     e                 max_var_v,max_var_a,max_param,param,
     e                 varint_v,varint_a,
     e                 varint_v_cero,varint_a_cero,max_def,
     e                 deftot_cero,deftot,esf_cero,esf,
     s                 deriv_g_esf)
c
c******************************************************************
c        LLAMADA A LA SUBRUTINA  QUE CALCULA LA
c        DERIVADA DE g CON RESPECTO A LAS VI CON FUENCION DE FLUENCIA
c******************************************************************
c        
        call cal_deriv_g_intv(max_prop,prop,max_opcion,opcion,
     e                 max_jprop,jprop,n_prop,max_int_a,max_int_v,
     e                 max_var_v,max_var_a,max_param,param,
     e                 varint_v,varint_a,
     e                 varint_v_cero,varint_a_cero,max_def,
     e                 deftot_cero,deftot,esf_cero,esf,
     s                 deriv_g_intv)
c
c******************************************************************
c        LLAMADA A LA SUBRUTINA  QUE CALCULA LA
c        DERIVADA DE g CON RESPECTO A LAS VI SIN FUENCION DE FLUENCIA
c******************************************************************
c        
        call cal_deriv_g_inta(max_prop,prop,max_opcion,opcion,
     e                 max_jprop,jprop,n_prop,max_int_a,
     e                 max_var_v,max_var_a,max_param,param,
     e                 varint_v,varint_a,
     e                 varint_v_cero,varint_a_cero,max_def,
     e                 deftot_cero,deftot,esf_cero,esf,
     s                 deriv_g_inta)
c
c******************************************************************
c                ENSAMBLAJE DE HIPERMATRIZ
c******************************************************************
c
        nfi=jprop(4)
        nv=jprop(5)
        na=jprop(6)
c
        do i=1,nfi
                do j=1,nfi
                        hiper(i,j)=deriv_r_esf(i,j)
                end do
c
                do j=nfi+1,nfi+nv
                        hiper(i,j)=deriv_r_intv(i,j-nfi)
                end do
c
                do j=nfi+nv+1,nfi+nv+na
                        hiper(i,j)=deriv_r_inta(i,j-nfi-nv)
                end do
        end do
c
        do i=nfi+1,nfi+nv
                do j=1,nfi
                        if(activas(i-nfi))then
                                hiper(i,j)=deriv_f_esf(i-nfi,j)
                        else
                                hiper(i,j)=0.
                        end if
                end do
c
                do j=nfi+1,nfi+nv
                        if(activas(i-nfi))then
                                hiper(i,j)=
     1                                deriv_f_intv(i-nfi,j-nfi)
                        else 
                            if(i.eq.j)then
                                hiper(i,j)=1.
                        else
                                     hiper(i,j)=0.
                      end if
                        end if
                end do
c
                do j=nfi+nv+1,nfi+nv+na
                        if(activas(i-nfi))then
                     hiper(i,j)=deriv_f_inta(i-nfi,j-nfi-nv)
                        else            
                                   hiper(i,j)=0.
                        end if
                end do
        end do
c
        do i=1+nfi+nv,nfi+nv+na
                do j=1,nfi
                        hiper(i,j)=deriv_g_esf(i-nfi-nv,j)
                end do
c
                do j=nfi+1,nfi+nv
                           hiper(i,j)=deriv_g_intv(i-nfi-nv,j-nfi)
                end do
c
                do j=nfi+nv+1,nfi+nv+na
                        hiper(i,j)=deriv_g_inta(i-nfi-nv,j-nfi-nv)
                end do
        end do       
c
        return
        end
c
c
c******************************************************************
c
c       SUBRUTINA  CAL_DERIV_R_ESF
c       CALCULO DE LAS DERIVADAS DE {R} CON RESPECTO A {M}
c       DONDE {R}={f-fp}-[F(D+)]{M+}-[F(D-)]{M-}
c
c
c       J. FLOREZ LOPEZ
c       A. RAMIREZ
c       (Nov./94)
c
c
c******************************************************************
c
        subroutine cal_deriv_r_esf(max_prop,prop,max_opcion,opcion,
     e         max_jprop,jprop,n_prop,max_esf,
     e         max_var_v,max_var_a,max_param,param,varint_v_cero,
     e         varint_a_cero,max_def,
     e         deftot_cero,deftot,esf_cero,
     e         varint_v,varint_a,esf,max_hiper,
     s         deriv_r_esf)
c
       IMPLICIT NONE
c******************************************************************
c               DECLARACION DE VARIABLES NUCLEO
c******************************************************************
c
c       prop= matriz de propiedades del elemento
c       param= parametros de la integracion numerica
c       opcion= opciones del calculo (dinamico o estatico,
c               peque~os o grandes desplazamientos, etc)
c       opcion(1)= "static" o "dynamic"
c       opcion(2)= "small" o "large"
c       dtiempo= incremento de tiempo en el paso
c       varint_a_cero= variables internas ai al principio del paso
c       varint_a= variables internas ai al final del paso
c       varint_v_cero= variables internas vi al principio del paso
c       varint_v= variables internas vi al final del paso
c       n_prop= numero de propiedades del elemento
c       max_prop= tama~o de la matriz de propiedades
c       max_var= tama~o de la matriz de variables internas
c       max_opcion= tama~o de la matriz de opciones
c       max_param= tama~o de la matriz de parametros
c       max_def= tama~o de la matriz de deformaciones
c       esf=esfuerzos al final del paso
c       deftot=deformaciones totales al final del paso
c       esf_cero= esfuerzos al principio del paso
c       deftot_cero= deformaciones totales al principio del paso
c       deftot_cero= deformaciones totales al principio del paso
c       activas= matriz que indica el estado de las variables internas
c
c*********************************************************************
c
        integer max_prop,max_opcion,max_param
        integer max_def,max_hiper,max_esf,max_jprop,jprop(max_jprop)
c
        Real(kind=8):: prop(max_prop),param(max_param)
c
        character(len=8):: opcion(max_opcion)
c
        integer max_var_a,max_var_v,n_prop
c
        Real(kind=8):: varint_a(max_var_a),varint_v(max_var_v),
     2               varint_a_cero(max_var_a),varint_v_cero(max_var_v)
c
        Real(kind=8):: deftot_cero(max_def),esf_cero(max_esf),
     1          deftot(max_def),esf(max_esf)
c
        Real(kind=8):: deriv_r_esf(max_esf,max_esf)
c
c*********************************************************************
c               DECLARACIONES LOCALES
c*********************************************************************
c
c       fpi= rotacion plastica en el nudo i
c       fpj= rotacion plastica en el nudo j
c       di_pos= da~o en el nudo i
c       di_neg
c       dj_pos= da~o en el nudo j
c       dj_neg
c       ei= modulo de elasticidad* momento de inercia=prop(4)
c       ea= modulo de elasticidad* area=prop(5)
c   
c       long= longitud del elemento
c       mi= momento en el extremo i
c       mj= momento en el extremo j
c       n= fuerza axial
c       fi= rotacion total en el extremo i
c       fj= rotacion total en el extremo j
c       fa= alargamiento de la cuerda
c       f11,f12,f21,f22,f33= elementos de la matriz de flexibilidad
c
c*********************************************************************
c
        Real(kind=8):: fpi,fpj,ei,ea,long,mi,mj,n,fi,fj,fa
        Real(kind=8):: di_pos,dj_pos,di_neg,dj_neg
c
c*********************************************************************
c               TRADUCCION NUCLEO --> MOD1
c*********************************************************************
c
        fi=deftot(1)
        fj=deftot(2)
        fa=deftot(3)
c
        mi=esf(1)
        mj=esf(2)
        n=esf(3)
c
        fpi=varint_v(1)
        fpj=varint_v(2)
        di_pos=varint_v(3)
        dj_pos=varint_v(4)
        di_neg=varint_v(5)
        dj_neg=varint_v(6)
c
        long=prop(1)
        ei=prop(4)
        ea=prop(5)
c
c*********************************************************************
c               CALCULO DE LAS DERIVADAS
c*********************************************************************
c       
        if(mi.ge.0.)then
                deriv_r_esf(1,1)=-long/(3.*ei*(1.-di_pos))
        else
                deriv_r_esf(1,1)=-long/(3.*ei*(1.-di_neg))
        end if
c
        deriv_r_esf(1,2)=long/(6.*ei)
        deriv_r_esf(1,3)=0.

        deriv_r_esf(2,1)=long/(6.*ei)
 
c
        if(mj.ge.0.)then
                deriv_r_esf(2,2)=-long/(3.*ei*(1.-dj_pos))
        else
                deriv_r_esf(2,2)=-long/(3.*ei*(1.-dj_neg))
        end if
c
        deriv_r_esf(2,3)=0.
c
        deriv_r_esf(3,1)=0.
        deriv_r_esf(3,2)=0.
c
        deriv_r_esf(3,3)=-long/ea
c
        return
        end
c
c******************************************************************
c
c        SUBRUTINA  CAL_DERIV_R_INTV
c        CALCULO DE LAS DERIVADAS DE {R} CON RESPECTO A 
c        LAS VARIABLES INTERNAS CON FUNCION DE FLUENCIA
c        DONDE {R}={f-fp}-[F(D+)]{M+}-[F(D-)]{M-}
c
c        ESCRITA POR J. FLOREZ LOPEZ y MARÕA EUGENIA MARANTE
c        (Mayo 2002)
c     Modificada por: Nayive Jaramillo S. (Enero 2002)
c
c******************************************************************
c
        subroutine cal_deriv_r_intv(max_prop,prop,max_opcion,opcion,
     e            max_jprop,jprop,n_prop,max_esf,max_int_v,
     e            max_var_v,max_var_a,max_param,param,varint_v_cero,
     e            varint_a_cero,max_def,
     e            deftot_cero,deftot,esf_cero,
     e            varint_v,varint_a,esf,max_hiper,
     s            deriv_r_intv)
c
       IMPLICIT NONE
c******************************************************************
c                DECLARACION DE VARIABLES NUCLEO
c******************************************************************
c
c        prop= matriz de propiedades del elemento
c        param= parametros de la integracion numerica
c        param(1)=alfa
c        param(2)=beta
c        param(3)=gamma
c        param(4)=mﬂximo nînero de iteraciones
c        param(5)=mﬂximo nîmero de pasos de integracin local
c       param(6)=alfa limite
c        param(7:6+nt)=ceros numericos para {R},{T},{X}
c        opcion= opciones del calculo (dinamico o estatico,
c                peque±os o grandes desplazamientos, etc)
c        opcion(1)= "static" o "dynamic"
c        opcion(2)= "small" o "large"
c        dtiempo= incremento de tiempo en el paso
c        varint_v_cero= variables internas con funcin de fluencia al principio del paso
c        varint_a_cero= variables internas sin funcin de fluencia al principio del paso
c        varint_v= variables internas con funcin de fluencia al final del paso
c        varint_a= variables internas sin funcin de fluencia al final del paso
c        nv= nîmero de variables internas con funcin de fluencia
c        na= nîmero de variables internas sin funcin de fluencia
c        ny= nîmero de fuerzas termodinﬂmicas
c        n_prop= nîmero de propiedades del elemento
c        nfi= dimensin de la matriz de deformaciones y de esfuerzos
c        max_prop= tama±o de la matriz de propiedades
c        max_var_v= tama±o de la matriz de var.int. con funcin de fluencia
c        max_var_a= tama±o de la matriz de var.int. sin funcin de fluencia 
c        max_opcion= tama±o de la matriz de opciones
c        max_param= tama±o de la matriz de parametros
c        max_def= tama±o de la matriz de deformaciones
c        esf=esfuerzos al final del paso
c        deftot=deformaciones totales al final del paso
c        esf_cero= esfuerzos al principio del paso
c        deftot_cero= deformaciones totales al principio del paso
c        fterm= Fuerzas termodinamicas asociadas a las variables 
c                    internas al final del paso
c        activas= matriz que indica el estado de todas las variables
c                  internas
c        max_act= tama±o de la matriz activas
c        rtg=({R},{T},{G})
c
c******************************************************************
c
        integer max_prop,max_var_v,max_var_a,max_opcion,max_param
c
        integer max_def,max_hiper,max_esf,max_int_v,max_jprop
c
        Real(kind=8):: prop(max_prop),param(max_param)
c
        character(len=8):: opcion(max_opcion)
c
        integer n_prop,jprop(max_jprop)
c
        Real(kind=8):: varint_v(max_var_v),varint_a(max_var_a),
     2       varint_v_cero(max_var_v),varint_a_cero(max_var_a)
c
        Real(kind=8):: deftot_cero(max_def),esf_cero(max_esf),
     1             deftot(max_def),esf(max_esf)
c
        Real(kind=8):: deriv_r_intv(max_esf,max_int_v)
c
*********************************************************************
c               DECLARACIONES LOCALES
c*********************************************************************
c
c       fpi= rotacion plastica en el nudo i
c       fpj= rotacion plastica en el nudo j
c       di_pos= da~o en el nudo i
c       di_neg
c       dj_pos= da~o en el nudo j
c       dj_neg
c       ei= modulo de elasticidad* momento de inercia=prop(4)
c       ea= modulo de elasticidad* area=prop(5)
c       c_pos= coeficiente del criterio de plasticidad
c       c_neg=                                        
c       my_pos= coeficiente del criterio de plasticidad
c       my_neg=                                       
c       long= longitud del elemento
c       mi= momento en el extremo i
c       mj= momento en el extremo j
c       n= fuerza axial
c       fi= rotacion total en el extremo i
c       fj= rotacion total en el extremo j
c       fa= alargamiento de la cuerda
c       f11,f12,f21,f22,f33= elementos de la matriz de felxibilidad
c
c*********************************************************************
c
        Real(kind=8):: fpi,fpj,ei,ea,long,mi,mj,n,fi,fj,fa
        Real(kind=8):: di_pos,dj_pos,di_neg,dj_neg
        Real(kind=8):: pos_mi,neg_mi,pos_mj,neg_mj
c
c*********************************************************************
c               TRADUCCION VARIABLES
c*********************************************************************
c
        fi=deftot(1)
        fj=deftot(2)
        fa=deftot(3)
c
        mi=esf(1)
        mj=esf(2)
        n=esf(3)
c
        fpi=varint_v(1)
        fpj=varint_v(2)
c
        di_pos=varint_v(3)
        dj_pos=varint_v(4)
        di_neg=varint_v(5)
        dj_neg=varint_v(6)
c
        long=prop(1)
        ei=prop(4)
        ea=prop(5)
c
c*********************************************************************
c               CALCULO DE LAS PARTES POSITIVA Y NEGATIVA DE LOS
C               MOMENTOS
c*********************************************************************
c
        if(mi.ge.0)then
                pos_mi=mi
                neg_mi=0.
        else
                pos_mi=0.
                neg_mi=mi
        end if
c
        if(mj.ge.0.)then
                pos_mj=mj
                neg_mj=0.
        else
                pos_mj=0.
                neg_mj=mj
        end if
c
c******************************************************************
c                CALCULO DE LAS DERIVADAS
c******************************************************************
c
        deriv_r_intv(1,1)=-1.
        deriv_r_intv(1,2)=0.
        deriv_r_intv(1,3)=-long*pos_mi/(3.*ei*(1.-di_pos)*(1.-di_pos))
        deriv_r_intv(1,4)=0.
        deriv_r_intv(1,5)=-long*neg_mi/(3.*ei*(1.-di_neg)*(1.-di_neg))
        deriv_r_intv(1,6)=0.
c
        deriv_r_intv(2,1)=0.
        deriv_r_intv(2,2)=-1.
        deriv_r_intv(2,3)=0.
        deriv_r_intv(2,4)=-long*pos_mj/(3.*ei*(1.-dj_pos)*(1.-dj_pos))
        deriv_r_intv(2,5)=0.
        deriv_r_intv(2,6)=-long*neg_mj/(3.*ei*(1.-dj_neg)*(1.-dj_neg))
c
        deriv_r_intv(3,1)=0.
        deriv_r_intv(3,2)=0.
        deriv_r_intv(3,3)=0.
        deriv_r_intv(3,4)=0.
        deriv_r_intv(3,5)=0.
        deriv_r_intv(3,6)=0.
c
        return
        end
c
c
c******************************************************************
c
c        SUBRUTINA  CAL_DERIV_R_INTA
c        CALCULO DE LAS DERIVADAS DE {R} CON RESPECTO A 
c        LAS VARIABLES INTERNAS SIN FUNCION DE FLUENCIA
c        DONDE {R}={f-fp}-[F(D+)]{M+}-[F(D-)]{M-}
c
c        ESCRITA POR J. FLOREZ LOPEZ y MARÕA EUGENIA MARANTE
c        (Mayo 2001)
c     Modificada por : Nayive Jaramillo S. (Enero 2002)
c******************************************************************
c
        subroutine cal_deriv_r_inta(max_prop,prop,max_opcion,opcion,
     e            max_jprop,jprop,n_prop,max_esf,max_int_a,
     e            max_var_v,max_var_a,max_param,param,varint_v_cero,
     e            varint_a_cero,max_def,
     e            deftot_cero,deftot,esf_cero,
     e            varint_v,varint_a,esf,max_hiper,
     s            deriv_r_inta)
c
       IMPLICIT NONE
c******************************************************************
c                DECLARACION DE VARIABLES NUCLEO
c******************************************************************
c
c        prop= matriz de propiedades del elemento
c        param= parametros de la integracion numerica
c        param(1)=alfa
c        param(2)=beta
c        param(3)=gamma
c        param(4)=mﬂximo nînero de iteraciones
c        param(5)=mﬂximo nîmero de pasos de integracin local
c       param(6)=alfa limite
c        param(7:6+nt)=ceros numericos para {R},{T},{X}
c        opcion= opciones del calculo (dinamico o estatico,
c                peque±os o grandes desplazamientos, etc)
c        opcion(1)= "static" o "dynamic"
c        opcion(2)= "small" o "large"
c        dtiempo= incremento de tiempo en el paso
c        varint_v_cero= variables internas con funcin de fluencia al principio del paso
c        varint_a_cero= variables internas sin funcin de fluencia al principio del paso
c        varint_v= variables internas con funcin de fluencia al final del paso
c        varint_a= variables internas sin funcin de fluencia al final del paso
c        nv= nîmero de variables internas con funcin de fluencia
c        na= nîmero de variables internas sin funcin de fluencia
c        ny= nîmero de fuerzas termodinﬂmicas
c        n_prop= nîmero de propiedades del elemento
c        nfi= dimensin de la matriz de deformaciones y de esfuerzos
c        max_prop= tama±o de la matriz de propiedades
c        max_var_v= tama±o de la matriz de var.int. con funcin de fluencia
c        max_var_a= tama±o de la matriz de var.int. sin funcin de fluencia 
c        max_opcion= tama±o de la matriz de opciones
c        max_param= tama±o de la matriz de parametros
c        max_def= tama±o de la matriz de deformaciones
c        esf=esfuerzos al final del paso
c        deftot=deformaciones totales al final del paso
c        esf_cero= esfuerzos al principio del paso
c        deftot_cero= deformaciones totales al principio del paso
c        fterm= Fuerzas termodinamicas asociadas a las variables 
c                    internas al final del paso
c        activas= matriz que indica el estado de todas las variables
c                  internas
c        max_act= tama±o de la matriz activas
c        rtg=({R},{T},{G})
c
c******************************************************************
c
        integer max_prop,max_var_v,max_var_a,max_opcion,max_param
c
        integer max_def,max_hiper,max_esf,max_int_a,max_jprop
c
        Real(kind=8):: prop(max_prop),param(max_param)
c
        character(len=8):: opcion(max_opcion)
c
        integer n_prop,jprop(max_jprop)
c
        Real(kind=8):: varint_v(max_var_v),varint_a(max_var_a),
     2           varint_v_cero(max_var_v),varint_a_cero(max_var_a)
c
        Real(kind=8):: deftot_cero(max_def),esf_cero(max_def),
     1                deftot(max_def),esf(max_def)
c
        Real(kind=8):: deriv_r_inta(max_esf,max_int_a)
c
c*********************************************************************
c               CALCULO DE LAS DERIVADAS
c*********************************************************************
c
        deriv_r_inta(1,1)=0.
        deriv_r_inta(1,2)=0.
c
        deriv_r_inta(2,1)=0.
        deriv_r_inta(2,2)=0.
c
        deriv_r_inta(3,1)=0.
        deriv_r_inta(3,2)=0.
C
        return
        end
c
c*********************************************************************
c
c        SUBRUTINA  CAL_DERIV_F_ESF
c        CALCULO DE LAS DERIVADAS DE {F} CON RESPECTO A {M} NUMERICAMENTE
c
c        DONDE {varint_v}=(fpi,fpij)
c              {f}=(f1,f2)
c
c        f1= funcin inelﬂstica asociada a fpi
c        f2= funcin inelﬂstica asociada a fpj
c
c        f_pos=  M - (1-d+)*alfa+*c+*fp - (1 - d+)*((1 - alfa+)*c+*p + My+)
c
c        f_neg= -M + (1-d-)*alfa-*c-*fp) - (1 - d-)*((1 - alfa- )*c-*p + My-)
c
c
c        ESCRITA POR Maria Eugenia Marante y Julio Florez Lopez (Mayo 2002)
c******************************************************************
c
        subroutine cal_deriv_f_esf(max_prop,prop,max_opcion,opcion,
     e              max_jprop,jprop,n_prop,max_esf,max_int_v,
     e              max_var_v,max_var_a,max_param,param,varint_v_cero,
     e              varint_a_cero,max_def,
     e              deftot_cero,deftot,esf_cero,
     e              varint_v,varint_a,esf,max_hiper,
     s              deriv_f_esf,conv)
c
       IMPLICIT NONE
c******************************************************************
c                DECLARACION DE VARIABLES NUCLEO
c******************************************************************
c
c        prop= matriz de propiedades del elemento
c        param= parametros de la integracion numerica
c        param(1)=alfa
c        param(2)=beta
c        param(3)=gamma
c        param(4)=mﬂximo nînero de iteraciones
c        param(5)=mﬂximo nîmero de pasos de integracin local
c        param(7:6+nt)=ceros numericos para {R},{T},{X}
c        opcion= opciones del calculo (dinamico o estatico,
c                peque±os o grandes desplazamientos, etc)
c        opcion(1)= "static" o "dynamic"
c        opcion(2)= "small" o "large"
c        dtiempo= incremento de tiempo en el paso
c        varint_v_cero= variables internas con funcin de fluencia al principio del paso
c        varint_a_cero= variables internas sin funcin de fluencia al principio del paso
c        varint_v= variables internas con funcin de fluencia al final del paso
c        varint_a= variables internas sin funcin de fluencia al final del paso
c        nv= nîmero de variables internas con funcin de fluencia
c        na= nîmero de variables internas sin funcin de fluencia
c        ny= nîmero de fuerzas termodinﬂmicas
c        n_prop= nîmero de propiedades del elemento
c        nfi= dimensin de la matriz de deformaciones y de esfuerzos
c        max_prop= tama±o de la matriz de propiedades
c        max_var_v= tama±o de la matriz de var.int. con funcin de fluencia
c        max_var_a= tama±o de la matriz de var.int. sin funcin de fluencia 
c        max_opcion= tama±o de la matriz de opciones
c        max_param= tama±o de la matriz de parametros
c        max_def= tama±o de la matriz de deformaciones
c        esf=esfuerzos al final del paso
c        deftot=deformaciones totales al final del paso
c        esf_cero= esfuerzos al principio del paso
c        deftot_cero= deformaciones totales al principio del paso
c        fterm= Fuerzas termodinamicas asociadas a las variables 
c                    internas al final del paso
c        activas= matriz que indica el estado de todas las variables
c                  internas
c        max_act= tama±o de la matriz activas
c        rtg=({R},{T},{G})
c
c******************************************************************
c
        integer max_prop,max_var_v,max_var_a,max_opcion,max_param
c
        integer max_def,max_hiper,max_esf,max_int_v,max_jprop
c
        Real(kind=8):: prop(max_prop),param(max_param)
c
        character(len=8):: opcion(max_opcion)
c
        integer n_prop,jprop(max_jprop)
c
        Real(kind=8):: varint_v(max_var_v),varint_a(max_var_a),
     1           varint_v_cero(max_var_v),varint_a_cero(max_var_a)
c
        Real(kind=8):: deftot_cero(max_def),esf_cero(max_def),
     1                deftot(max_def),esf(max_def)
c
        Real(kind=8):: deriv_f_esf(max_int_v,max_esf)
        logical conv
c
c******************************************************************
c                DECLARACIONES LOCALES
c******************************************************************
c
c        fpi= rotacion plastica en el nudo i
c        fpj= rotacion plastica en el nudo j
c        di= da~o en el nudo i
c        dj= da~o en el nudo j
c        ei= modulo de elasticidad* momento de inercia=prop(10)
c        ea= modulo de elasticidad* area=prop(11)
c        c= coeficiente del criterio de plasticidad=prop(12)
c        my= coeficiente del criterio de plasticidad=prop(13)
c        gcr=coeficiente del criterio de plasticidad=prop(14)
c        q=coeficiente del criterio de plasticidad=prop(15)
c        long= longitud del elemento
c        mi= momento en el extremo i
c        mj= momento en el extremo j
c        n= fuerza axial
c        fi= rotacion total en el extremo i
c        fj= rotacion total en el extremo j
c        fa= alargamiento de la cuerda
c        f_plus= variacion de las funciones de fluencia
c        esf_plus= variacion de los esfuerzos
c        mp= momento plastico
c
c*********************************************************************
c
        integer i,j,max_f_plus,max_def_plus
          integer nfi,nv
        parameter (max_f_plus=6)
        parameter (max_def_plus=3)
        Real(kind=8):: f_plus(max_f_plus)
        Real(kind=8):: esf_plus(max_def_plus),f(max_f_plus)
        Real(kind=8):: my
c
        Real(kind=8):: valor_pos,valor_neg
          Real(kind=8):: fpi,fpj,mi,mj
        Real(kind=8):: c_pos_i,my_pos_i
        Real(kind=8):: c_pos_j,my_pos_j
        Real(kind=8):: c_neg_i,my_neg_i
        Real(kind=8):: c_neg_j,my_neg_j
        Real(kind=8):: di_pos,dj_pos,di_neg,dj_neg
        Real(kind=8):: alfa_pos_i,alfa_neg_i
        Real(kind=8):: alfa_pos_j,alfa_neg_j
c
c******************************************************************
c                TRADUCCION VARIABLES
c******************************************************************
c                
          my=prop(8)
        nfi=jprop(4)
        nv=jprop(5)
c
        mi=esf(1)
        mj=esf(2)
c
        fpi=varint_v(1)
        fpj=varint_v(2)
        di_pos=varint_v(3)
        dj_pos=varint_v(4)
        di_neg=varint_v(5)
        dj_neg=varint_v(6)
c
        c_pos_i=prop(6)
        c_neg_i=prop(14)
        c_pos_j=prop(10)
        c_neg_j=prop(18)
c
        alfa_pos_i=prop(55)
        alfa_neg_i=prop(55)
        alfa_pos_j=prop(55)
        alfa_neg_j=prop(55)
c
c
c******************************************************************
c                CALCULO DE LAS DERIVADAS POR DIFERENCIAS FINITAS
c******************************************************************
c
        call cal_f(max_prop,prop,max_opcion,opcion,
     e          max_var_v,max_var_a,n_prop,
     e          max_param,param,varint_v_cero,varint_a_cero,
     e          max_def,max_jprop,jprop,
     e          deftot_cero,deftot,esf_cero,
     e          varint_v,varint_a,esf,max_int_v,
     s          f,conv)
c
c
         if (.not. conv) then
        return
         endif

        do j=1,nfi-1        
              do i=1,nfi
                   esf_plus(i)=esf(i)
              end do
c
              esf_plus(j)=esf_plus(j)+
     1                         0.00001*my*0.77
c
                call cal_f(max_prop,prop,max_opcion,opcion,
     e          max_var_v,max_var_a,n_prop,
     e          max_param,param,varint_v_cero,varint_a_cero,
     e          max_def,max_jprop,jprop,
     e          deftot_cero,deftot,esf_cero,
     e          varint_v,varint_a,esf_plus,max_int_v,
     s          f_plus,conv)
c
               if (.not. conv) then
                return
        endif

              do i=1,nv
                   deriv_f_esf(i,j)=(f_plus(i)-f(i))/
     1                      (esf_plus(j)-esf(j))
              end do
        end do
c
        do i=1,nv
               deriv_f_esf(i,3)=0.
        end do
c
c
c
c******************************************************************
c        correccion de valores conocidos
c******************************************************************
c
c
          valor_pos=mi/(1.-di_pos)/(1.-di_neg)-alfa_pos_i*c_pos_i*fpi
        valor_neg=-mi/(1.-di_pos)/(1.-di_neg)+alfa_neg_i*c_neg_i*fpi
c
          if(valor_pos.ge.valor_neg)then
                    deriv_f_esf(1,1)=1.
          else
                    deriv_f_esf(1,1)=-1.
          end if
c
          deriv_f_esf(1,2)=0.
          deriv_f_esf(1,3)=0.
c
c************************************************************************
c                
          deriv_f_esf(2,1)=0.
c
          valor_pos=mj/(1.-dj_pos)/(1.-dj_neg)-alfa_pos_j*c_pos_j*fpj
        valor_neg=-mj/(1.-dj_pos)/(1.-dj_neg)+alfa_neg_j*c_neg_j*fpj
c
          if(valor_pos.ge.valor_neg)then
                    deriv_f_esf(2,2)=1.
          else
                    deriv_f_esf(2,2)=-1.
          end if
c
          deriv_f_esf(2,3)=0.
c
c************************************************************************
c                
c
        return
        end
c
c******************************************************************
c
c        SUBRUTINA  CAL_DERIV_G_ESF
c        CALCULO DE LAS DERIVADAS DE {G} CON RESPECTO A {M}
c
c
c
c        Escrita por Maria Eugenia Marante y Julio Florez Lopez
c                        (mayo2002)
c
c******************************************************************
c
        subroutine cal_deriv_g_esf(max_prop,prop,max_opcion,opcion,
     e                   max_jprop,jprop,n_prop,max_int_a,max_esf,
     e                   max_var_v,max_var_a,max_param,param,
     e                   varint_v,varint_a,
     e                   varint_v_cero,varint_a_cero,max_def,
     e                   deftot_cero,deftot,esf_cero,esf,
     s                   deriv_g_esf)
c
       IMPLICIT NONE
c******************************************************************
c                DECLARACION DE VARIABLES NUCLEO
c******************************************************************
c
        integer max_prop,max_var_v,max_var_a,max_opcion,max_param
        integer max_def,max_esf,max_int_a,max_jprop
        integer n_prop,jprop(max_jprop)
c
        character(len=8):: opcion(max_opcion)
c
        Real(kind=8):: prop(max_prop),param(max_param)
        Real(kind=8):: varint_v(max_var_v),varint_v_cero(max_var_v)
        Real(kind=8):: varint_a(max_var_a),varint_a_cero(max_var_a)
        Real(kind=8):: deftot_cero(max_def),esf_cero(max_esf)
             Real(kind=8):: deftot(max_def),esf(max_esf)
        Real(kind=8):: deriv_g_esf(max_int_a,max_esf)
c
c******************************************************************
c
        deriv_g_esf(1,1)=0.
        deriv_g_esf(1,2)=0.
        deriv_g_esf(1,3)=0.
        deriv_g_esf(2,1)=0.
        deriv_g_esf(2,2)=0.
        deriv_g_esf(2,3)=0.
c
        return
        end
c
c
c******************************************************************
c
c        SUBRUTINA  CAL_DERIV_F_INTV
c        CALCULO DE LAS DERIVADAS DE {F} CON RESPECTO A LAS
c        VARIABLES INTERNAS CON FUNCION DE FLUENCIA
C
c        DONDE {varint_v}=(fpi,fpj,pi,pj)
c              {f}=(f1,f2,f3,f4)
c
c        f1= funcin inelﬂstica asociada a fpi
c        f2= funcin inelﬂstica asociada a fpj
c        f3= funcin inelﬂstica asociada a pi
c        f4= funcin inelﬂstica asociada a pj
c
c        ESCRITA POR MARÕA EUGENIA MARANTE y JULIO FLOREZ LOPEZ 
c        (mayo 2001)
c     
c        Nayive Jaramillo (Enero 2002)
c******************************************************************
c
        subroutine cal_deriv_f_intv(max_prop,prop,max_opcion,opcion,
     e        max_jprop,jprop,n_prop,max_int_v,max_int_a,max_var_v,
     e        max_var_a,max_param,param,varint_v_cero,varint_a_cero,
     e        max_def,deftot_cero,deftot,esf_cero,
     e        varint_v,varint_a,esf,max_hiper,
     s        deriv_f_intv,conv)
c
       IMPLICIT NONE
c******************************************************************
c                DECLARACION DE VARIABLES NUCLEO
c******************************************************************
c
c        prop= matriz de propiedades del elemento
c        param= parametros de la integracion numerica
c        param(1)=alfa
c        param(2)=beta
c        param(3)=gamma
c        param(4)=mﬂximo nînero de iteraciones
c        param(5)=mﬂximo nîmero de pasos de integracin local
c        param(7:6+nt)=ceros numericos para {R},{T},{X}
c        opcion= opciones del calculo (dinamico o estatico,
c                peque±os o grandes desplazamientos, etc)
c        opcion(1)= "static" o "dynamic"
c        opcion(2)= "small" o "large"
c        dtiempo= incremento de tiempo en el paso
c        varint_v_cero= variables internas con funcin de fluencia al principio del paso
c        varint_a_cero= variables internas sin funcin de fluencia al principio del paso
c        varint_v= variables internas con funcin de fluencia al final del paso
c        varint_a= variables internas sin funcin de fluencia al final del paso
c        nv= nîmero de variables internas con funcin de fluencia
c        na= nîmero de variables internas sin funcin de fluencia
c        ny= nîmero de fuerzas termodinﬂmicas
c        n_prop= nîmero de propiedades del elemento
c        nfi= dimensin de la matriz de deformaciones y de esfuerzos
c        max_prop= tama±o de la matriz de propiedades
c        max_var_v= tama±o de la matriz de var.int. con funcin de fluencia
c        max_var_a= tama±o de la matriz de var.int. sin funcin de fluencia 
c        max_opcion= tama±o de la matriz de opciones
c        max_param= tama±o de la matriz de parametros
c        max_def= tama±o de la matriz de deformaciones
c        esf=esfuerzos al final del paso
c        deftot=deformaciones totales al final del paso
c        esf_cero= esfuerzos al principio del paso
c        deftot_cero= deformaciones totales al principio del paso
c        fterm= Fuerzas termodinamicas asociadas a las variables 
c                    internas al final del paso
c        activas= matriz que indica el estado de todas las variables
c                  internas
c        max_act= tama±o de la matriz activas
c        rtg=({R},{T},{G})
c
c******************************************************************
c
        integer max_prop,max_var_v,max_var_a,max_opcion,max_param
c
        integer max_def,max_hiper,max_int_v,max_int_a,max_jprop
c
        Real(kind=8):: prop(max_prop),param(max_param)
c
        character(len=8):: opcion(max_opcion)
c
        integer n_prop,jprop(max_jprop)
c
        Real(kind=8):: varint_v(max_var_v),varint_a(max_var_a),
     1           varint_v_cero(max_var_v),varint_a_cero(max_var_a)
c
        Real(kind=8):: deftot_cero(max_def),esf_cero(max_def),
     1                deftot(max_def),esf(max_def)
c
        Real(kind=8):: deriv_f_intv(max_int_v,max_int_v)
        logical conv

c*********************************************************************
c                DECLARACIONES LOCALES
c*********************************************************************
c        fpi= rotacion plastica en el nudo i
c        fpj= rotacion plastica en el nudo j
c        di= da~o en el nudo i
c        dj= da~o en el nudo j
c        ei= modulo de elasticidad* momento de inercia=prop(10)
c        ea= modulo de elasticidad* area=prop(11)
c        c= coeficiente del criterio de plasticidad=prop(12)
c        my= coeficiente del criterio de plasticidad=prop(13)
c        gcr=coeficiente del criterio de plasticidad=prop(14)
c        q=coeficiente del criterio de plasticidad=prop(15)
c        long= longitud del elemento
c        mi= momento en el extremo i
c        mj= momento en el extremo j
c        n= fuerza axial
c        fi= rotacion total en el extremo i
c        fj= rotacion total en el extremo j
c        fa= alargamiento de la cuerda
c     mp= momento plastico de la seccion a un nivel de carga axial
c        f_plus= variacion de las funciones de fluencia
c        varint_v_plus= variacion de las variables internas vi
c
c*********************************************************************
c
        integer i,j,max_int_v_plus
        parameter (max_int_v_plus=6)
        Real(kind=8):: f_plus(max_int_v_plus)
        Real(kind=8):: varint_v_plus(max_int_v_plus),f(max_int_v_plus)
          integer nv
c
c*********************************************************************
c                TRADUCCION VARIABLES
c*********************************************************************
         nv=jprop(5)
c
c
c******************************************************************
c                CALCULO DE LAS DERIVADAS POR DIFERENCIAS FINITAS
c******************************************************************
c
         call cal_f(max_prop,prop,max_opcion,opcion,
     e            max_var_v,max_var_a,n_prop,
     e            max_param,param,varint_v_cero,varint_a_cero,
     e            max_def,max_jprop,jprop,
     e                    deftot_cero,deftot,esf_cero,   
     e            varint_v,varint_a,esf,max_int_v,
     s            f,conv)
c
        if (.not. conv) then
        return
        endif
       do j=1,nv 
c       
            do i=1,nv
c
               varint_v_plus(i)=varint_v(i)
c
            end do
c
            if(varint_v(j).eq.0.)then
                  varint_v_plus(j)=0.0001
            else
                  varint_v_plus(j)=varint_v(j)+0.0001*varint_v(j)
            end if
c
c
          call cal_f(max_prop,prop,max_opcion,opcion,
     e                 max_var_v,max_var_a,n_prop,
     e                 max_param,param,varint_v_cero,varint_a_cero,
     e                 max_def,max_jprop,jprop,
     e                       deftot_cero,deftot,esf_cero,   
     e                 varint_v_plus,varint_a,esf,max_int_v,
     s                 f_plus,conv)
c
        if (.not. conv) then
                return
        endif
                  do i=1,nv
c
                      deriv_f_intv(i,j)=((f_plus(i)-f(i))/
     1               (varint_v_plus(j)-varint_v(j)))
            end do
        end do                                        
c
c
c******************************************************************
C        correccion de valores conocidos
c******************************************************************
c
c        deriv_f_intv(3,1)=0.
          deriv_f_intv(3,2)=0.
c
        deriv_f_intv(4,1)=0.
          deriv_f_intv(4,2)=0.
c
        deriv_f_intv(5,1)=0.
          deriv_f_intv(5,2)=0.
c
        deriv_f_intv(6,1)=0.
          deriv_f_intv(6,2)=0.
c
        return
        end
c
c******************************************************************
c
c        SUBRUTINA  CAL_DERIV_F_INTA
c        CALCULO DE LAS DERIVADAS DE {F} CON RESPECTO A LAS
c        VARIABLES INTERNAS SIN FUNCION DE FLUENCIA
c
c        DONDE {varint_a}=(di+,dj+,di-,dj-)
c              {f}=(f1,f2,f3,f4)
c
c        f1= funcin inelﬂstica asociada a fpi
c        f2= funcin inelﬂstica asociada a fpj
c 
c
c        ESCRITA POR MARÕA EUGENIA MARANTE y JULIO FLOREZ LOPEZ 
c        (mayo 2001)
c
c        Modificada por: Nayive Jaramillo (ENERO 2002)
c******************************************************************
c
        subroutine cal_deriv_f_inta(max_prop,prop,max_opcion,opcion,
     e        max_jprop,jprop,n_prop,max_int_v,max_int_a,max_var_v,
     e        max_var_a,max_param,param,varint_v_cero,varint_a_cero,
     e        max_def,deftot_cero,deftot,esf_cero,
     e        varint_v,varint_a,esf,max_hiper,
     s        deriv_f_inta,conv)
c
       IMPLICIT NONE
c******************************************************************
c                DECLARACION DE VARIABLES NUCLEO
c******************************************************************
c
c        prop= matriz de propiedades del elemento
c        param= parametros de la integracion numerica
c        param(1)=alfa
c        param(2)=beta
c        param(3)=gamma
c        param(4)=mﬂximo nînero de iteraciones
c        param(5)=mﬂximo nîmero de pasos de integracin local
c        param(7:6+nt)=ceros numericos para {R},{T},{X}
c        opcion= opciones del calculo (dinamico o estatico,
c                peque±os o grandes desplazamientos, etc)
c        opcion(1)= "static" o "dynamic"
c        opcion(2)= "small" o "large"
c        dtiempo= incremento de tiempo en el paso
c        varint_v_cero= variables internas con funcin de fluencia al principio del paso
c        varint_a_cero= variables internas sin funcin de fluencia al principio del paso
c        varint_v= variables internas con funcin de fluencia al final del paso
c        varint_a= variables internas sin funcin de fluencia al final del paso
c        nv= nîmero de variables internas con funcin de fluencia
c        na= nîmero de variables internas sin funcin de fluencia
c        ny= nîmero de fuerzas termodinﬂmicas
c        n_prop= nîmero de propiedades del elemento
c        nfi= dimensin de la matriz de deformaciones y de esfuerzos
c        max_prop= tama±o de la matriz de propiedades
c        max_var_v= tama±o de la matriz de var.int. con funcin de fluencia
c        max_var_a= tama±o de la matriz de var.int. sin funcin de fluencia 
c        max_opcion= tama±o de la matriz de opciones
c        max_param= tama±o de la matriz de parametros
c        max_def= tama±o de la matriz de deformaciones
c        esf=esfuerzos al final del paso
c        deftot=deformaciones totales al final del paso
c        esf_cero= esfuerzos al principio del paso
c        deftot_cero= deformaciones totales al principio del paso
c        fterm= Fuerzas termodinamicas asociadas a las variables 
c                    internas al final del paso
c        activas= matriz que indica el estado de todas las variables
c                  internas
c        max_act= tama±o de la matriz activas
c        rtg=({R},{T},{G})
c
c******************************************************************
c
        integer max_prop,max_var_v,max_var_a,max_opcion,max_param
c
        integer max_def,max_hiper,max_int_v,max_int_a,max_jprop
c
        Real(kind=8):: prop(max_prop),param(max_param)
c
        character(len=8):: opcion(max_opcion)
c
        integer n_prop,jprop(max_jprop)
c
        Real(kind=8):: varint_v(max_var_v),varint_a(max_var_a),
     1           varint_v_cero(max_var_v),varint_a_cero(max_var_a)
c
        Real(kind=8):: deftot_cero(max_def),esf_cero(max_def),
     1                deftot(max_def),esf(max_def)
c
        Real(kind=8):: deriv_f_inta(max_int_v,max_int_a)
        logical conv
c
c*********************************************************************
c                DECLARACIONES LOCALES
c*********************************************************************
c
c
c        fpi= rotacion plastica en el nudo i
c        fpj= rotacion plastica en el nudo j
c        di= da~o en el nudo i
c        dj= da~o en el nudo j
c        ei= modulo de elasticidad* momento de inercia=prop(10)
c        ea= modulo de elasticidad* area=prop(11)
c        c= coeficiente del criterio de plasticidad=prop(12)
c        my= coeficiente del criterio de plasticidad=prop(13)
c        gcr=coeficiente del criterio de plasticidad=prop(14)
c        q=coeficiente del criterio de plasticidad=prop(15)
c        long= longitud del elemento
c        mi= momento en el extremo i
c        mj= momento en el extremo j
c        n= fuerza axial
c        fi= rotacion total en el extremo i
c        fj= rotacion total en el extremo j
c        fa= alargamiento de la cuerda
c       mp= momento plastico de la seccion a un nivel de carga axial
c        h_plus= variacion de las funciones de fluencia, deslizamiento y da~os
c        var_int_plus= variacion de las variables internas
c
c*********************************************************************
c
        integer i,j,max_int_v_plus,max_int_a_plus
        parameter (max_int_v_plus=6)
        parameter (max_int_a_plus=2)
        Real(kind=8):: f_plus(max_int_v_plus)
        Real(kind=8):: varint_a_plus(max_int_a_plus),f(max_int_v_plus)
          integer na,nv
c
c*********************************************************************
c                TRADUCCION NUCLEO --> MDC
c*********************************************************************
c
       na=jprop(6)
         nv=jprop(5)
c
c*********************************************************************
c
c******************************************************************
c                CALCULO DE LAS DERIVADAS POR DIFERENCIAS FINITAS
c******************************************************************
c
c
         call cal_f(max_prop,prop,max_opcion,opcion,
     e            max_var_v,max_var_a,n_prop,
     e            max_param,param,varint_v_cero,varint_a_cero,
     e            max_def,max_jprop,jprop,
     e                deftot_cero,deftot,esf_cero,   
     e            varint_v,varint_a,esf,max_int_v,
     s            f,conv)
c
        if (.not. conv) then
                return
        endif

          do j=1,na
           do i=1,na
c
                varint_a_plus(i)=varint_a(i)
c
            end do
c
            if(varint_a(j).eq.0.)then
                  varint_a_plus(j)=0.0001
            else
                  varint_a_plus(j)=varint_a(j)+0.0001*varint_a(j)
            end if
c
                  call cal_f(max_prop,prop,max_opcion,opcion,
     e                 max_var_v,max_var_a,n_prop,
     e                 max_param,param,varint_v_cero,varint_a_cero,
     e                 max_def,max_jprop,jprop,
     e                       deftot_cero,deftot,esf_cero,   
     e                 varint_v,varint_a_plus,esf,max_int_v,
     s                 f_plus,conv)
c
c
            if (.not. conv) then
                return
            endif        
                  do i=1,nv
c
                deriv_f_inta(i,j)=((f_plus(i)-f(i))/
     1               (varint_a_plus(j)-varint_a(j)))
            end do
        end do                                        
c
c
c******************************************************************
C        correccion de valores conocidos
c******************************************************************
c
c          deriv_f_inta(3,1)=0.
          deriv_f_inta(3,2)=0.
          deriv_f_inta(4,1)=0.
          deriv_f_inta(4,2)=0.
          deriv_f_inta(5,1)=0.
          deriv_f_inta(5,2)=0.
          deriv_f_inta(6,1)=0.
          deriv_f_inta(6,2)=0.
c
        return
        end
c
c******************************************************************
c
c        SUBRUTINA  CAL_DERIV_G_INTV
c        CALCULO DE LAS DERIVADAS DE {G} CON RESPECTO A {INTV}
c
c        donde  {g}=(g1,g2)
c                   {varint_v}=(fpi,fpj,di+,dj+,di-,dj-)
c        
c
c        ESCRITA POR JULIO FL‘REZ L‘PEZ y MARÕA EUGENIA MARANTE
c        (Mayo 2001)
c     Modificada por:
c     NAYIVE JARAMILLO
c     (ENERO 2002)
c
c******************************************************************
c
        subroutine cal_deriv_g_intv(max_prop,prop,max_opcion,opcion,
     e                   max_jprop,jprop,n_prop,max_int_a,max_int_v,
     e                   max_var_v,max_var_a,max_param,param,
     e                   varint_v,varint_a,
     e                   varint_v_cero,varint_a_cero,max_def,
     e                   deftot_cero,deftot,esf_cero,esf,
     s                   deriv_g_intv)
c
       IMPLICIT NONE
c******************************************************************
c                DECLARACION DE VARIABLES NUCLEO
c******************************************************************
c
        integer max_prop,max_var_v,max_var_a,max_opcion,max_param
        integer max_def,max_int_a,max_int_v,max_jprop
        integer n_prop,jprop(max_jprop)
c
        character(len=8):: opcion(max_opcion)
c
        Real(kind=8):: prop(max_prop),param(max_param)
        Real(kind=8):: varint_v(max_var_v),varint_v_cero(max_var_v)
        Real(kind=8):: varint_a(max_var_a),varint_a_cero(max_var_a)
        Real(kind=8):: deftot_cero(max_def),esf_cero(max_def)
             Real(kind=8):: deftot(max_def),esf(max_def)
        Real(kind=8):: deriv_g_intv(max_int_a,max_int_v)
        Real(kind=8):: fpi,fpj,pi_cero,pj_cero
c
c*********************************************************************
c
        fpi=varint_v(1)
        fpj=varint_v(2)
        pi_cero=varint_a_cero(1)
        pj_cero=varint_a_cero(2)
c
c******************************************************************
c
        if(dabs(fpi).ge.pi_cero)then
                if(fpi.ge.0.)then
                        deriv_g_intv(1,1)=-1.
                else
                        deriv_g_intv(1,1)=1.
                end if
        else
                deriv_g_intv(1,1)=0.
        end if
c
        deriv_g_intv(1,2)=0.
        deriv_g_intv(1,3)=0.
        deriv_g_intv(1,4)=0.
        deriv_g_intv(1,5)=0.
        deriv_g_intv(1,6)=0.
c
c
        deriv_g_intv(2,1)=0.
c
        if(dabs(fpj).ge.pj_cero)then
                if(fpj.ge.0.)then
                        deriv_g_intv(2,2)=-1.
                else
                        deriv_g_intv(2,2)=1.
                end if
        else
                deriv_g_intv(2,2)=0.
        end if                        
c
        deriv_g_intv(2,3)=0.
        deriv_g_intv(2,4)=0.
        deriv_g_intv(2,5)=0.
        deriv_g_intv(2,6)=0.
c
      return
      end
c
c******************************************************************
c
c        SUBRUTINA  CAL_DERIV_G_INTA
c        CALCULO DE LAS DERIVADAS DE {G} CON RESPECTO A {INTA}
c
c        donde  {g}=(g1,g2)
c
c                 {varint_a}=(pi,pj)
c        
c
c        ESCRITA POR JULIO FL‘REZ L‘PEZ y MARÕA EUGENIA MARANTE
c        (Mayo 2001)
c     Modificada por:
c     NAYIVE JARAMILLO
c     (ENERO 2002)
c******************************************************************
c
        subroutine cal_deriv_g_inta(max_prop,prop,max_opcion,opcion,
     e                        max_jprop,jprop,n_prop,max_int_a,
     e                        max_var_v,max_var_a,max_param,param,
     e                        varint_v,varint_a,
     e                        varint_v_cero,varint_a_cero,max_def,
     e                        deftot_cero,deftot,esf_cero,esf,
     s                        deriv_g_inta)
c
       IMPLICIT NONE
c******************************************************************
c                DECLARACION DE VARIABLES NUCLEO
c******************************************************************
c
        integer max_prop,max_var_v,max_var_a,max_opcion,max_param
        integer max_def,max_int_a,max_jprop
        integer n_prop,jprop(max_jprop)
c
        character(len=8):: opcion(max_opcion)
c
        Real(kind=8):: prop(max_prop),param(max_param)
        Real(kind=8):: varint_v(max_var_v),varint_v_cero(max_var_v)
        Real(kind=8):: varint_a(max_var_a),varint_a_cero(max_var_a)
        Real(kind=8):: deftot_cero(max_def),esf_cero(max_def)
             Real(kind=8):: deftot(max_def),esf(max_def)
        Real(kind=8):: deriv_g_inta(max_int_a,max_int_a)
c
c*********************************************************************
c
          deriv_g_inta(1,1)=1.
          deriv_g_inta(1,2)=0.
c
          deriv_g_inta(2,1)=0.
          deriv_g_inta(2,2)=1.
c
        return
        end
c
c*********************************************************************
c*********************************************************************
c********                                                  ***********
c********     FACTORIZACION Y SOLUCION DEL SISTEMA DE      ***********
c********              ECUACIONES NO LINEAL                ***********
c********                                                  ***********
c********  UTILIZANDO EL PROCESO DE ELIMINACION DE GAUSS   ***********
c********                Y PIVOTE MAXIMO                   ***********
c********                      (Ax=b)                      ***********
c*********************************************************************
c
c      ESCRITA EN EL 2001 POR:
c      ING. DENIS AVON
c      ING. RICARDO PICON
c
c*********************************************************************
c
c
        subroutine rfssedp(a,ind,n)
c
       IMPLICIT NONE
c*********************************************************************
c
c               DECLARACIONES LOCALES
c
c*********************************************************************
c
c        a: es la matriz asimâtrica (A) que se va a factorizar.
c
c        ind: es el vector (b) independiente y tambiân serﬂ el vector
c             solucin (x) al final del proceso.
c
c        n: es el orden de la matriz cuadrada asimâtrica.
c
c*********************************************************************
c
        integer n,I,L,ip,J,K,bandera,t
        Real(kind=8):: RMAX,a(n,n),aux1,aux2,ind(n),rm
c
c
c********************************************************************
c
c      COMIENZA EL PROCESO DE ELIMINACION HACIA ADELANTE
c
c********************************************************************
c
c
        do K=1,(n-1)
c
c            BUSQUEDA DEL PIVOTE MAXIMO
c
            RMAX=0.
            bandera=0
            do J=K,n
                if(abs(a(J,K)).gt.RMAX)then
                     ip=J
                     RMAX=abs(a(J,K))
                     bandera=1
                end if
            end do
c
c*************************************************
c         PRIMER CHEQUEO DE LA SOLUCION
c           SISTEMA INDETERMINADO
c*************************************************
c
            if(bandera.eq.0)then
                  do t=1,n
                     ind(t)=0.
                  end do
                  return
            end if
c
c*************************************************
c         CAMBIO DE FILAS
c*************************************************
c
c
            if(ip.gt.K)then
                 do L=K,n
                     aux1=a(K,L)
                     a(K,L)=a(ip,L)
                     a(ip,L)=aux1
                 end do
                 aux2=ind(K)
                 ind(K)=ind(ip)
                 ind(ip)=aux2
            end if
c
c***************************************************
c                ELIMINACION DE INCOGNITAS
c***************************************************
c
c
            do I=K+1,n
                 rm=a(I,K)/a(K,K)
                 do J=K+1,n
                     a(I,J)=a(I,J)-rm*a(K,J)
                 end do
                 ind(I)=ind(I)-rm*ind(K)
            end do
        end do
c
c**************************************************
c         SEGUNDO CHEQUEO DE LA SOLUCION
c           SISTEMA INDETERMINADO
c**************************************************
c
        if(abs(a(n,n)).lt.0.000001)then
c
                  do t=1,n
                     ind(t)=0.
                  end do
                  return
        end if
c
c**************************************************
c       CALCULO DE LA ULTIMA INCOGNITA
c**************************************************
c
        ind(n)=ind(n)/a(n,n)
c
c******************************************************
c     COMIENZA EL PROCESO DE SUSTITUCION HACIA ATRAS
c******************************************************
c
c
        do L=1,n-1
            K=n-L
            do J=K+1,n
                ind(K)=ind(K)-a(K,J)*ind(J)
            end do
            ind(K)=ind(K)/a(K,K)
        end do
c
c
        return
        end
c
c
c*****************************************************************
c        SUBRUTINA QUE ACTUALIZA LOS ESFUERZOS Y LAS VARIABLES
c        INTERNAS
c*****************************************************************
c
        subroutine actualizacion(max_prop,prop,max_jprop,jprop,
     e                 n_prop,max_var_v,max_var_a,max_param,param,
     e                 max_def,max_hiper,rtg,
     s                 esf,varint_v,varint_a)
c
       IMPLICIT NONE
c******************************************************************
c        DECLARACION DE VARIABLES MDC
c******************************************************************
c
c        prop= matriz de propiedades del elemento
c        param= parametros de la integracion numerica
c        param(1)=alfa
c        param(2)=beta
c        param(3)=gamma
c        param(4)=maximo numero de iteraciones
c        param(5)=maximo numero de pasos de integracin local
c        param(7:6+nfi+nv+na)=ceros numâricos de {R},{T},{G}
c        nv= nîmero de variables internas con funcin de fluencia
c        na= nîmero de variables internas sin funcin de fluencia
c        nfi= dimensin de la matriz de deformaciones y de esfuerzos
c        n_prop= numero de propiedades del elemento
c        max_prop= tama~o de la matriz de propiedades
c        max_var_v= tama~o de la matriz de variables internas con funcin de fluencia
c        max_param= tama~o de la matriz de parametros
c        max_def= tama~o de la matriz de deformaciones
c        esf=esfuerzos al final del paso
c
c******************************************************************
c
        integer max_prop,max_var_v,max_var_a,max_def,max_param,max_jprop
c
        Real(kind=8):: prop(max_prop),param(max_param)
c
        integer n_prop,max_hiper,jprop(max_jprop)
c
        Real(kind=8):: varint_v(max_var_v),varint_a(max_var_a)
c
        Real(kind=8):: esf(max_def),rtg(max_hiper)
c
c******************************************************************
c                DECLARACIONES LOCALES
c******************************************************************
c
        integer i,nfi,nv,na
c
c******************************************************************
c
        nfi=jprop(4)
        nv=jprop(5)
        na=jprop(6)
c
c
        do i=1,nfi
                esf(i)=esf(i)-rtg(i)
        end do
c
        do i=1,nv
                varint_v(i)=varint_v(i)-rtg(i+nfi)
        end do
c
        do i=1,na
                varint_a(i)=varint_a(i)-rtg(i+nfi+nv)
        end do
c
        return 
        end
c
c
c******************************************************************
c
c	SUBRUTINA  CONVER
c        DETERMINA LA CONVERGENCIA DEL SIGUIENTE SISTEMA DE
c        ECUACIONES NO LINEALES:
c
c        R(esfuerzos,varint_v,varint_a,deftot)=0
c        T= (Vr - Vr-1)j fj(esfuerzos,varint_v,varint_a,deftot)=0 donde f<=0
c        G(esfuerzos,varint_v,varint_a,deftot)=0
c
c        ....
c        DONDE
c
c        R representa a la ley de estado
c        T representa la ley de evolucin de las variables internas Vi
c        G representa la Ley de Dan~o que define las variables internas Ai 
c        f es la funcion de fluencia de la variable interna Vj
c
c
c        ESCRITA POR JULIO FLOREZ LOPEZ y MARÕA EUGENIA MARANTE
c        (Abril 2001)
c
c     Modificada por : Nayive Jaramillo S.(ENERO 2002)
c******************************************************************
c
        subroutine conver(rtg,param,max_param,max_hiper,max_jprop,
     e                        jprop,max_def,cambio,
     s                        convergencia)
c
       IMPLICIT NONE
c******************************************************************
c                DECLARACION DE VARIABLES ESFUERZOS
c******************************************************************
c
c        rtv= ({R},{T},{G})= funciones que deben igualarse a cero
c        max_hiper= tama~o maximo de rtg
c        param= parametros de la integracion numerica
c        param(1)=alfa
c        param(2)=beta
c        param(3)=gamma
c        param(4)=maximo numero de iteraciones
c        param(5)=maximo numero de pasos de integracin =PNEWDT
c        param(6)=alfa limite
c        param(7:6+nfi+nv+na)=ceros numâricos de {R},{T},{G}
c        max_param= tama~o de la matriz de parametros
c        convergencia= variable logica que indica la convergencia o no
c                          del proceso de resolucion del sistema de ecuaciones
c        nv= numero de variables internas con funcin de fluencia
c        na= nîmero de variables internas sin funcin de fluencia
c        nfi= dimensin de la matriz de deformaciones y de esfuerzos
c        max_def= tama~o de la matriz de deformaciones
c        cambio= variable logica que indica un cambio de activas a pasivas
c                    o viceversa en la subrutina verif no en residual3d
c
c******************************************************************
c
        integer max_param,max_hiper,max_jprop,max_def
c
        integer jprop(max_jprop)
c
        Real(kind=8):: rtg(max_hiper),param(max_param)
c
        logical convergencia,cambio
c
c******************************************************************
c                DECLARACIONES LOCALES
c******************************************************************
c
c        tol= ceros numericos para {R,T,G}
c
c******************************************************************
c
        integer i,nfi,nv,na,max_tol
        parameter(max_tol=11)
        Real(kind=8):: tol(max_tol)
c
c
        nfi=jprop(4)
        nv=jprop(5)
        na=jprop(6)
c
c
        do i=1,nfi
                tol(i)=param(6+i)
        end do
c
        do i=nfi+1,nfi+nv
                tol(i)=param(6+i)
        end do
c
        do i=nfi+nv+1,nfi+nv+na
                tol(i)=param(6+i)
        end do
c
        convergencia=.true.
c
        if (cambio)then
                convergencia=.false.
                cambio=.false.
        end if
c
        i=1
        do while(convergencia.and.i.le.nfi+nv+na)
                if(dabs(rtg(i)).gt.tol(i))then
                        convergencia=.false.
                end if
                i=i+1
        end do
c
        return
        end
c
c******************************************************************
c
c        SUBRUTINA VERIF
c                
c        ESCRITA POR J. FLOREZ LOPEZ y MARÕA EUGENIA MARANTE
c        (Abril 2001)
c     Modificada por : Nayive Jaramillo (Enero 2002)
c******************************************************************
c
        subroutine verif(max_prop,prop,max_opcion,opcion,
     e           max_jprop,jprop,n_prop,dtiempo,
     e           max_var_v,max_var_a,max_param,param,varint_v_cero,
     e           varint_a_cero,max_def,deftot_cero,deftot,esf_cero,
     e           varint_v,varint_a,esf,max_act,max_f,f,
     s           activas,cambio)
c
       IMPLICIT NONE
c******************************************************************
c                DECLARACION DE VARIABLES NUCLEO
c******************************************************************
c
c        prop= matriz de propiedades del elemento
c        param= parametros de la integracion numerica
c        param(1)=alfa
c        param(2)=beta
c        param(3)=gamma
c        param(4)=mﬂximo nînero de iteraciones
c        param(5)=mﬂximo nîmero de pasos de integracin=PNEWDT
c        param(6)=alafa limite
c        param(7:6+nt)=ceros numericos para {R},{T},{G}
c        opcion= opciones del calculo (dinamico o estatico,
c                peque±os o grandes desplazamientos, etc)
c        opcion(1)= "static" o "dynamic"
c        opcion(2)= "small" o "large"
c        dtiempo= incremento de tiempo en el paso
c        varint_v_cero= variables internas con funcin de fluencia al principio del paso
c        varint_a_cero= variables internas sin funcin de fluencia al principio del paso
c        varint_v= variables internas con funcin de fluencia al final del paso
c        varint_a= variables internas sin funcin de fluencia al final del paso
c        nv= nîmero de variables internas con funcin de fluencia
c        na= nîmero de variables internas sin funcin de fluencia
c        ny= nîmero de fuerzas termodinﬂmicas
c        n_prop= nîmero de propiedades del elemento
c        nfi= dimensin de la matriz de deformaciones y de esfuerzos
c        max_prop= tama±o de la matriz de propiedades
c        max_var_v= tama±o de la matriz de var.int. con funcin de fluencia
c        max_var_a= tama±o de la matriz de var.int. sin funcin de fluencia 
c        max_opcion= tama±o de la matriz de opciones
c        max_param= tama±o de la matriz de parametros
c        max_def= tama±o de la matriz de deformaciones
c        esf=esfuerzos al final del paso
c        deftot=deformaciones totales al final del paso
c        esf_cero= esfuerzos al principio del paso
c        deftot_cero= deformaciones totales al principio del paso
c        fterm= Fuerzas termodinamicas asociadas a las variables 
c                    internas al final del paso
c        activas= matriz que indica el estado de todas las variables
c                  internas
c        max_act= tama±o de la matriz activas
c        rtg=({R},{T},{G})
c        cambio= variable logica que indica un cambio de activa a pasiva 
c                de alguna variable interna
c
c******************************************************************
c
        integer max_prop,max_var_v,max_var_a,max_opcion,max_param
c
        integer max_jprop,max_def,max_act,max_f
c
        Real(kind=8):: prop(max_prop),param(max_param),
     2                dtiempo
c
        character(len=8):: opcion(max_opcion)
c
        integer n_prop,jprop(max_jprop)
c
        Real(kind=8):: varint_v(max_var_v),varint_a(max_var_a),
     2           varint_v_cero(max_var_v),varint_a_cero(max_var_a)
c
        Real(kind=8):: deftot_cero(max_def),esf_cero(max_def),
     1                deftot(max_def),esf(max_def),f(max_f)
c
        logical activas(max_act),cambio
c
c*********************************************************************
c               DECLARACIONES LOCALES
c*********************************************************************
c
      integer i,nfi,nv,na,nt
      Real(kind=8):: tol(max_var_a+max_var_v+max_def)
c
        Real(kind=8):: mi,mj
        Real(kind=8):: fpi,fpj,fpi_cero,fpj_cero
        Real(kind=8):: di_pos,dj_pos,di_neg,dj_neg
        Real(kind=8):: di_pos_cero,dj_pos_cero,di_neg_cero,dj_neg_cero
        Real(kind=8):: alfa_pos_i,alfa_neg_i,alfa_pos_j,alfa_neg_j
        Real(kind=8):: c_pos_i,c_neg_i,c_pos_j,c_neg_j
        Real(kind=8):: valor_pos_i,valor_neg_i
        Real(kind=8):: valor_pos_j,valor_neg_j
c
c*********************************************************************
c
      mi=esf(1)
      mj=esf(2)
c
      fpi=varint_v(1)
      fpj=varint_v(2)
c
      fpi_cero=varint_v_cero(1)
      fpj_cero=varint_v_cero(2)
c
      di_pos=varint_v(3)
      dj_pos=varint_v(4)
      di_neg=varint_v(5)
      dj_neg=varint_v(6)
c
      di_pos_cero=varint_v_cero(3)
      dj_pos_cero=varint_v_cero(4)
      di_neg_cero=varint_v_cero(5)
      dj_neg_cero=varint_v_cero(6)
c
      alfa_pos_i=prop(55)
      alfa_neg_i=prop(55)
      alfa_pos_j=prop(55)
      alfa_neg_j=prop(55)
c
      c_pos_i=prop(6)
      c_neg_i=prop(14)
      c_pos_j=prop(10)
      c_neg_j=prop(18)
c                                
c*********************************************************************
c               VERIFICACION DE LAS DEFORMACIONES PLASTICAS
c*********************************************************************
c
       valor_pos_i=mi/(1.-di_pos)/(1.-di_neg)-alfa_pos_i*c_pos_i*fpi
c
       valor_neg_i=-mi/(1.-di_pos)/(1.-di_neg)+alfa_neg_i*c_neg_i*fpi
c
       valor_pos_j=mj/(1.-dj_pos)/(1.-dj_neg)-alfa_pos_j*c_pos_j*fpj
c
       valor_neg_j=-mj/(1.-dj_pos)/(1.-dj_neg)+alfa_neg_j*c_neg_j*fpj
c
c
       if(activas(1).and.valor_pos_i.gt.valor_neg_i.and.
     1   (fpi-fpi_cero).lt.0.)then
c
              activas(1) = .false.
                cambio = .true.
         end if
c
      if(activas(1).and.valor_pos_i.lt.valor_neg_i.and.
     1   (fpi-fpi_cero).gt.0.)then
c
              activas(1) = .false.
                cambio = .true.
         end if
c        

       if(activas(2).and.valor_pos_j.gt.valor_neg_j.and.
     1   (fpj-fpj_cero).lt.0.)then
c
              activas(2) = .false.
                cambio = .true.
         end if
c
       if(activas(2).and.valor_pos_j.lt.valor_neg_j.and.
     1   (fpj-fpj_cero).gt.0.)then
c
              activas(2) = .false.
                cambio = .true.
         end if
c
        if (activas(3).and.(di_pos-di_pos_cero).lt.0.) then
                activas(3) = .false.
                cambio = .true.
        end if
c
        if (activas(4).and.(dj_pos-dj_pos_cero).lt.0.) then
                activas(4) = .false.
                cambio = .true.
        end if
c
        if (activas(5).and.(di_neg-di_neg_cero).lt.0.) then
                activas(5) = .false.
                cambio = .true.
        end if
c
        if (activas(6).and.(dj_neg-dj_neg_cero).lt.0.) then
                activas(6) = .false.
                cambio = .true.
        end if
c
      return
      end
c
c******************************************************************
c
c        SUBRUTINA  ACTIV
c     ESTA SUBRUTINA VERIFICA SI LAS VARIABLES INTERNAS
c     CUMPLEN CON LAS CONDICIONES TERMODINAMICAS.
c
c        ESCRITA POR Marça Eugenia Marante y Julio Flrez Lpez
c        Marzo 2001
c     Modificada por: NAYIVE JARAMILLO (Enero 2002)
c******************************************************************
c
        subroutine activ(max_prop,prop,
     e            max_jprop,jprop,n_prop,
     e            max_var_v,max_param,param,max_def,
     e            max_act,max_f,f,
     s            activas,cambio)
c
       IMPLICIT NONE
c******************************************************************
c        DECLARACION DE VARIABLES MDC1
c******************************************************************
c
c        prop= matriz de propiedades del elemento
c        param= parametros de la integracion numerica
c        param(1)=alfa
c        param(2)=beta
c        param(3)=gamma
c        param(4)=maximo numero de iteraciones
c        param(5)=maximo numero de pasos de integracin =PNEWDT
c        param(7:6+nfi+nv+na)=ceros numâricos de {R},{T},{G}
c        nv= nîmero de variables internas con funcin de fluencia
c        nfi= dimensin de la matriz de deformaciones y de esfuerzos
c        n_prop= numero de propiedades del elemento
c        max_prop= tama~o de la matriz de propiedades
c        max_var_v= tama~o de la matriz de variables internas con funcin de fluencia
c        max_param= tama~o de la matriz de parametros
c        max_def= tama~o de la matriz de deformaciones
c        activas= matriz que indica el estado de las variables internas
c
c******************************************************************
c
        integer max_prop,max_jprop,max_var_v,max_param
c
        integer max_def,max_act,max_f
c
        Real(kind=8):: prop(max_prop),param(max_param)
c
        integer n_prop,jprop(max_jprop),nfi
c
        logical activas(max_act),cambio
c
        Real(kind=8):: f(max_f)
c
c******************************************************************
c                DECLARACIONES LOCALES
c******************************************************************
c
        integer i
c
        integer max_fi
c
        parameter (max_fi=11)
c        
        Real(kind=8):: tol(max_fi)
c
        integer nv,nf
c
c******************************************************************
c
        nfi=jprop(4)
        nv=jprop(5)
c
c******************************************************************
c        CEROS NUMERICOS PARA LAS FUNCIONES INELASTICAS
c******************************************************************
c
        do i=nfi+1,nfi+nv
                tol(i)=param(i+6)
        end do
c
        do i=1,nv
                if(.not.activas(i).and.f(i).gt.tol(i+nfi))then
                        activas(i)=.true.
                        cambio=.true.
                end if
        end do
c
        return
        end
c
c******************************************************************
c
c        SUBRUTINA  CAL_FTERM
c        CALCULO DE LAS FUERZAS TERMODINAMICAS ASOCIADAS
c        A LAS VARIABLES INTERNAS
c
c          ESCRITA POR J. FLOREZ LOPEZ y MARÕA EUGENIA MARANTE
c          (Abril 2001) 
c       Modificada por:NAYIVE JARAMILLO
c       (ENERO 2001)
c
c******************************************************************
c
        subroutine cal_fterm(max_prop,prop,max_opcion,opcion,
     e           max_fterm,max_jprop,jprop,n_prop,dtiempo,
     e           max_var_v,max_var_a,max_param,param,varint_v_cero,
     e           max_def,deftot_cero,deftot,esf_cero,
     e           varint_v,varint_a,esf,
     s           fterm)
c
       IMPLICIT NONE
c******************************************************************
c                DECLARACION DE VARIABLES NUCLEO
c******************************************************************
c
c        prop= matriz de propiedades del elemento
c        param= parametros de la integracion numerica
c        param(1)=alfa
c        param(2)=beta
c        param(3)=gamma
c        param(4)=mﬂximo nînero de iteraciones
c        param(5)=mﬂximo nîmero de pasos de integracin =PNEWDT
c        param(7:6+nt)=ceros numericos para {R},{T},{G}
c        opcion= opciones del calculo (dinamico o estatico,
c                peque±os o grandes desplazamientos, etc)
c        opcion(1)= "static" o "dynamic"
c        opcion(2)= "small" o "large"
c        dtiempo= incremento de tiempo en el paso
c        varint_v_cero= variables internas con funcin de fluencia al principio del paso
c        varint_a_cero= variables internas sin funcin de fluencia al principio del paso
c        varint_v= variables internas con funcin de fluencia al final del paso
c        varint_a= variables internas sin funcin de fluencia al final del paso
c        nv= nîmero de variables internas con funcin de fluencia
c        na= nîmero de variables internas sin funcin de fluencia
c        ny= nîmero de fuerzas termodinﬂmicas
c        n_prop= nîmero de propiedades del elemento
c        nfi= dimensin de la matriz de deformaciones y de esfuerzos
c        max_prop= tama±o de la matriz de propiedades
c        max_var_v= tama±o de la matriz de var.int. con funcin de fluencia
c        max_var_a= tama±o de la matriz de var.int. sin funcin de fluencia 
c        max_opcion= tama±o de la matriz de opciones
c        max_param= tama±o de la matriz de parametros
c        max_def= tama±o de la matriz de deformaciones
c        esf=esfuerzos al final del paso
c        deftot=deformaciones totales al final del paso
c        esf_cero= esfuerzos al principio del paso
c        deftot_cero= deformaciones totales al principio del paso
c        fterm= Fuerzas termodinamicas asociadas a las variables 
c                    internas al final del paso
c        activas= matriz que indica el estado de todas las variables
c                  internas
c        max_act= tama±o de la matriz activas
c        rtg=({R},{T},{G})
c        cambio= variable logica que indica un cambio de activa a pasiva 
c                de alguna variable interna
c
c******************************************************************
c
        integer max_prop,max_var_v,max_var_a,max_opcion,max_param
        integer max_fterm
        integer max_def,max_jprop
        integer n_prop,jprop(max_jprop)
c
        character(len=8):: opcion(max_opcion)
c
        Real(kind=8):: prop(max_prop),param(max_param),dtiempo
        Real(kind=8):: varint_v(max_var_v),varint_a(max_var_a)
        Real(kind=8):: varint_v_cero(max_var_v)
        Real(kind=8):: deftot_cero(max_def),esf_cero(max_def)
        Real(kind=8):: deftot(max_def),esf(max_def)
        Real(kind=8):: fterm(max_fterm)
c
c****c*********************************************************************
c
        Real(kind=8):: ei,long,mi,mj
        Real(kind=8):: di_pos,dj_pos,di_neg,dj_neg
        Real(kind=8):: pos_mi,neg_mi,pos_mj,neg_mj
        Real(kind=8):: f0,Gi_pos,Gj_pos,Gi_neg,Gj_neg
        Real(kind=8):: n,fpi,fpj
c
c*********************************************************************
c
        ei=prop(4)
        long=prop(1)
c
        mi=esf(1)
        mj=esf(2)
        n=esf(3)
c
        fpi=varint_v(1)
        fpj=varint_v(2)
c
        di_pos=varint_v(3)
        dj_pos=varint_v(4)
        di_neg=varint_v(5)
        dj_neg=varint_v(6)
c
c*********************************************************************
c               CALCULO DE LAS PARTES POSITIVA Y NEGATIVA DE LOS
C               MOMENTOS
c*********************************************************************
c
        if(mi.ge.0)then
                pos_mi=mi
                neg_mi=0.
        else
                pos_mi=0.
                neg_mi=mi
        end if
c
        if(mj.ge.0.)then
                pos_mj=mj
                neg_mj=0.
        else
                pos_mj=0.
                neg_mj=mj
        end if
c
        f0=long/(3*ei)
        Gi_pos=pos_mi*pos_mi*f0/(2.*(1.-di_pos)*(1.-di_pos))
        Gj_pos=pos_mj*pos_mj*f0/(2.*(1.-dj_pos)*(1.-dj_pos))
c
        Gi_neg=neg_mi*neg_mi*f0/(2.*(1.-di_neg)*(1.-di_neg))
        Gj_neg=neg_mj*neg_mj*f0/(2.*(1.-dj_neg)*(1.-dj_neg))
c
        fterm(1)=Gi_pos
        fterm(2)=Gj_pos
        fterm(3)=Gi_neg
        fterm(4)=Gj_neg
c
        return
        end
c
c******************************************************************
c        Subrutina para la gestion del paso local
c        escrita por Denis Avon 10/05/2001
c        modificada por Maria Eugenia Marante  y Julio Florez-Lopez
c        08/10/2001
c     Modificada por: Nayive Jaramillo (ENERO 2002)
c*******************************************************************
c
        subroutine cambio_de_paso(max_def,max_deftot_int,max_var_v,
     e           max_var_a,deftot,deftot_int,varint_v,varint_a,esf,
     s           varint_v_cero,varint_a_cero,esf_cero,max_param,
     s           param,convergencia,deftot_cero,alpha)
c
       IMPLICIT NONE
c*******************************************************************
c        Declaraciones
c*******************************************************************
c
        integer max_deftot_int,max_def,max_var_v,max_var_a
        integer max_param,i
        Real(kind=8):: alpha,alpha_lim,deftot_int(max_deftot_int)
        Real(kind=8):: deftot_cero(max_def),deftot(max_def)
        Real(kind=8):: varint_v(max_var_v),varint_a(max_var_a)
        Real(kind=8):: varint_v_cero(max_var_v),varint_a_cero(max_var_a)
        Real(kind=8):: esf_cero(max_def),param(max_param),esf(max_def)
        logical convergencia
c
c**********************************************************************
c
      alpha_lim=1./10000.
c        
        if (.not.convergencia) then
            alpha=alpha/2.
            if (alpha.lt.alpha_lim) then             
c
                param(5)=0.25
                convergencia=.true.
                      alpha=1.
                return
c
            else
                return
            end if
        else if (alpha.eq.1.) then
            return
        else
            do i=1,max_deftot_int
                 deftot_cero(i)=deftot_int(i)
            end do
c
            do i=1,max_def
                 esf_cero(i)=esf(i)
            end do
c
            do i=1,max_var_v
                 varint_v_cero(i)=varint_v(i)
            end do        
c
            do i=1,max_var_a
                 varint_a_cero(i)=varint_a(i)
            end do        
c
            alpha=1.
            convergencia=.false.
        end if
c
        return
        end
c
c******************************************************************
c
c        SUBRUTINA  CAL_RESIDU
c        CALCULO DE LAS FUERZAS RESIDUALES AL FINAL DEL PASO
c
c        
c        VERSION 1.0
c        ESCRITA POR Marça Eugenia Marante y Julio Flrez Lpez
c        Abril 2001
c
c     Modificada por: Nayive Jaramillo (ENERO 2002)
c******************************************************************
c
        subroutine cal_residu(max_prop,prop,max_desp,desp,max_opcion,
     e                  opcion,param,max_param,esf_cero,
     e                  max_def,esf,acel,velo,
     e                  nq,coord,max_nud,max_coo,
     s                  residu)
c
       IMPLICIT NONE
c******************************************************************
c        DECLARACION DE VARIABLES MDC
c******************************************************************
c
c        prop= matriz de propiedades del elemento
c        desp= desplazamientos nodales al final del paso
c        velo= velocidades nodales al final del paso
c        acel= aceleraciones nodales al final del paso
c        param= parametros de la integracion numerica
c        opcion= opciones del calculo (dinamico o estatico,
c                peque~os o grandes desplazamientos, etc)
c        opcion(1)= "static" o "dynamic"
c        opcion(2)= "small" o "large"
c        dtiempo= incremento de tiempo en el paso
c        fuerza= fuerzas internas al final del paso
c        residu= fuerzas residuales al final del paso
c        max_prop= tama~o de la matriz de propiedades
c        max_opcion= tama~o de la matriz de opciones
c        max_param= tama~o de la matriz de parametros
c        max_desp= numero maximo de desplazamientos
c        coord= coordenadas de los nudos del elemento
c        max_nud= maximo numero de nudos del elemento=2
c        max_coo= maximo numero de coordenadas=2
c        max_def= tama~o de la matriz de deformaciones
c        esf=esfuerzos al final del paso
c        esf_cero=esfuerzos en el paso anterior
c        nq= nîmero de grados de libertad
c
c******************************************************************
c
        integer max_prop,max_desp,max_opcion,max_param
c
        integer nq,max_nud,max_coo,max_def,i
c
        Real(kind=8):: prop(max_prop),desp(max_desp),velo(max_desp),
     1         acel(max_desp)
c
        character(len=8):: opcion(max_opcion)
c
        Real(kind=8):: coord(max_coo,max_nud),param(max_param)
c
        Real(kind=8):: esf(max_def),residu(max_desp),
     1         esf_cero(max_def)
c
c******************************************************************
c        DECLARACIONES LOCALES
c******************************************************************
c
c        long= longitud del elemento
c       seno= seno del angulo con la horizontal
c       coseno= coseno del angulo con la horizontal
c        fuerza= matriz de fuerzas internas
c        fuerza_inerc= matriz de fuerzas de inercia
c        bt= traspuesta de la matriz de trasformacion
c        masa= matriz de masas
c        max_lib= max_desp
c        max_esf= max_def
c       alfa=parametro de integracion numerica para el caso dinamico
c        fuerza_cero=fuerzas internas en el paso anterior
c
c******************************************************************
c
        integer max_lib,max_esf
        parameter(max_lib=6)
        parameter(max_esf=3)
        Real(kind=8):: bt(max_lib,max_esf),b(max_esf,max_lib)
        Real(kind=8):: fuerza(max_lib),fuerza_inerc(max_lib)
        Real(kind=8):: masa(max_lib,max_lib),long,coseno,seno
        Real(kind=8):: fuerza_cero(max_lib),alfa
c
c******************************************************************
c                CALCULO DE LAS FUERZAS INTERNAS  
c******************************************************************
c
c        CALCULO DE LA TRASPUESTA DE LA MATRIZ DE
c        TRANSFORMACION
c
c
        long=prop(1)
        seno=prop(2)
        coseno=prop(3)
        call cal_b(long,seno,coseno,
     e                  max_def,max_lib,
     s                  b)
c
        call cal_bt(b,max_def,max_lib,
     s              bt)
c
c        CALCULO DE LAS FUERZAS INTERNAS
c
        call pro_2mat(max_desp,max_def,max_def,1,
     1                max_desp,1,
     1                max_desp,max_def,1,
     2                bt,esf,
     1                fuerza)
c
c
c******************************************************************
c        CALCULO DE LAS FUERZAS DE INERCIA 
c        TODAVIA NO SE CONSIDERAN LAS FUERZAS AMORTIGUAMIENTO
c******************************************************************
c
         if (opcion(1).eq."dynamic")then
c
              alfa=param(1)
c
c******************************************************************
c        CALCULO DE LA MATRIZ DE MASA DEL MIEMBRO
c******************************************************************
c
              call cal_masa(prop,max_prop,max_coo,max_nud,
     e                  max_desp,max_def,coord,
     s                  masa)
c
c******************************************************************
c        CALCULO DE LA MATRIZ DE FUERZAS DE INERCIA
c******************************************************************
c
              call pro_2mat(max_desp,max_desp,
     e                  max_desp,1,
     e                  max_desp,1,
     e                  max_desp,max_desp,1,
     e                  masa,acel,
     s                  fuerza_inerc)
c
c******************************************************************
c        CALCULO DE LAS FUERZAS DE AMORTIGUAMIENTO NUMERICO
c        ver manual ABAQUS-THEORY 3.4.1
c        A.RAMIREZ  8/9/94
c******************************************************************
c
             call pro_2mat(max_desp,max_def,max_def,1,
     e                 max_desp,1,max_desp,max_def,1,
     e                 bt,esf_cero,
     s                 fuerza_cero)
c
             do i=1,nq
                    residu(i)=fuerza_inerc(i)+
     1                            (alfa+1.)*fuerza(i)-
     2                             alfa*fuerza_cero(i)
             end do
c
        else
             if(opcion(1).eq."static".or.opcion(1).eq."ini_acel")then
                 do i=1,nq
                       residu(i)=fuerza(i)
                 end do
c
             else
                 stop
     1               "error en la matriz opcion"
             end if
c
        end if
c
        return
        end        
c        
c
c******************************************************************
c
c        SUBRUTINA  CAL_JACOB
c        CALCULO DEL JACOBIANO GLOBAL
c
c        ESCRITA POR JULIO FLOREZ LOPEZ y MARÕA EUGENIA MARANTE
c        (Abril 2001)
c     MODIFICADA POR: NAYIVE JARAMILLO (ENERO 2002)
c******************************************************************
c
        subroutine cal_jacob(max_prop,prop,max_desp,desp,max_opcion,
     e             max_jprop,jprop,opcion,nq,n_prop,max_var_v,
     e             max_var_a,max_def,max_coo,max_nud,esf,dtiempo,
     e             coord,deftot,varint_v,varint_a,
     e             param,max_param,max_act,activas,
     e             varint_v_cero,varint_a_cero,deftot_cero,esf_cero,
     s             jacob)
c
       IMPLICIT NONE
c******************************************************************
c        DECLARACION DE VARIABLES MDC
c******************************************************************
c
c        prop= matriz de propiedades del elemento
c        param= parametros de la integracion numerica
c        param(1)=alfa
c        param(2)=beta
c        param(3)=gamma
c        param(4)=maximo numero de iteraciones
c        param(5)=maximo numero de pasos de integracin=PNEWDT
c        param(7:6+nfi+nv+na)=ceros numâricos de {R},{T},{G}
c        opcion= opciones del calculo (dinamico o estatico,
c                peque~os o grandes desplazamientos, etc)
c        opcion(1)= "static" o "dynamic"
c        opcion(2)= "small" o "large"
c        dtiempo= incremento de tiempo en el paso
c        varint_v_cero= variables internas con funcin de fluencia al principio del paso
c        varint_v = variables internas con funcin de fluencia al final del paso
c        varint_a_cero= variables internas sin funcin de fluencia al principio del paso
c        varint_a = variables internas sin funcin de fluencia al final del paso
c        nq=  nîmero de grados de libertad
c        nv= nîmero de variables internas con funcin de fluencia
c        na= nîmero de variables internas sin funcin de fluencia
c        nfi= dimensin de la matriz de deformaciones y de esfuerzos
c        n_prop= numero de propiedades del elemento
c        max_prop= tama~o de la matriz de propiedades
c        max_var_v= tama~o de la matriz de variables internas con funcin de fluencia
c        max_var_a= tama~o de la matriz de variables internas sin funcin de fluencia
c        max_opcion= tama~o de la matriz de opciones
c        max_param= tama~o de la matriz de parametros
c        max_def= tama~o de la matriz de deformaciones
c        esf=esfuerzos al final del paso
c        deftot=deformaciones totales al final del paso
c        activas= matriz que indica el estado de las variables internas
c        hiper= hipermatriz para la resolucion del sistema de ecuaciones
c                locales por el metodo de Newton
c        max_hiper= tama~o maximo de hiper
c
c******************************************************************
c
        integer max_prop,max_var_v,max_var_a,max_opcion,max_param
c
        integer max_jprop,max_def,max_act,max_desp,max_coo,max_nud
c
        Real(kind=8):: prop(max_prop),param(max_param)
c
        character(len=8):: opcion(max_opcion)
c
        integer nq,n_prop,jprop(max_jprop)
c
        Real(kind=8):: varint_v(max_var_v),varint_a(max_var_a),
     e                 desp(max_desp)
c
        Real(kind=8):: deftot(max_def),esf(max_def),dtiempo
c
        Real(kind=8):: jacob(max_desp,max_desp),coord(max_coo,max_nud)
c
        Real(kind=8):: esf_cero(max_def),deftot_cero(max_def)
c
        Real(kind=8):: varint_v_cero(max_var_v),varint_a_cero(max_var_a)
c
        logical activas(max_act)
c
c******************************************************************
c                DECLARACIONES LOCALES
c******************************************************************
c
c        tang= matriz tangente local
c        max_esf=tama~o de la matriz de esfuerzos=max_def
c        b= matriz de transformacion
c        bt= traspuesta de la matriz de transformacion
c        masa= matriz de masas del elemento
c        max_for= max_desp
c        beta= parametro de la integracion numerica
c        deriva_u= derivadas de la aceleracion con respecto al desplazamiento
c        long= longitud del elemento
c        
c******************************************************************
c
        integer max_esf,max_for,i,j
        parameter(max_esf=3)
        parameter(max_for=6)
c
        Real(kind=8):: tang(max_esf,max_esf),b(max_esf,max_for)
        Real(kind=8):: bt(max_for,max_esf),masa(max_for,max_for)
        Real(kind=8):: deriva_u(max_for,max_for)
        Real(kind=8):: jacob_iner(max_for,max_for)
        Real(kind=8):: long,seno,coseno,beta
c
        Real(kind=8):: n
c
c******************************************************************
c                CALCULO DEL JACOBIANO LOCAL
c******************************************************************
c
        call cal_tang(max_prop,prop,max_jprop,jprop,max_opcion,
     e                        opcion,n_prop,max_var_v,max_var_a,
     e                        max_def,esf,varint_v_cero,varint_a_cero,
     e                        deftot,varint_v,varint_a,
     e                        esf_cero,deftot_cero,
     e                        param,max_param,max_act,activas,
     s                        tang)
c
c******************************************************************
c        CALCULO DE LA CONTRIBUCION LOCAL INTERNA AL
c        JACOBIANO GLOBAL
c******************************************************************
c        
c        CALCULO DE LA MATRIZ DE TRANSFORMACION Y DE SU 
c        TRASPUESTA 
c
        long=prop(1)
        seno=prop(2)
        coseno=prop(3)
c
        call cal_b(long,seno,coseno,
     e               max_def,max_desp,                      
     s                   b)
c
        call cal_bt(b,max_def,max_desp,
     s                        bt)
c
        call pro_3mat(max_desp,max_def,max_def,max_def,
     e                        max_def,max_desp,max_desp,max_desp,
     e                        max_desp,max_def,max_def,max_desp,
     e                        bt,tang,b,
     e                        jacob)
c
      if(opcion(2).eq."large")then
c
c----------------------------------------------------------------------
C
c       CONTRIBUCION DE LOS GRANDES DESPLAZAMIENTOS 
C       AL JACOBIANO  GLOBAL.
C       OPCION ESCRITA POR: ALEXIS LOPEZ INOJOSA
C       REALIZADA EL 17/2/94
c----------------------------------------------------------------------
c               
c
           n=esf(3)
c
           jacob(1,1)=jacob(1,1)+(n*seno**2)/long**2
           jacob(1,2)=jacob(1,2)-(n*seno*coseno)/long**2
           jacob(1,4)=jacob(1,4)-(n*seno**2)/long**2
           jacob(1,5)=jacob(1,5)+(n*seno*coseno)/long**2
           jacob(2,1)=jacob(2,1)-(n*seno*coseno)/long**2
           jacob(2,2)=jacob(2,2)+(n*coseno**2)/long**2
           jacob(2,4)=jacob(2,4)+(n*seno*coseno)/long**2
           jacob(2,5)=jacob(2,5)-(n*coseno**2)/long**2
           jacob(4,1)=jacob(4,1)-(n*seno**2)/long**2
           jacob(4,2)=jacob(4,2)+(n*seno*coseno)/long**2
           jacob(4,4)=jacob(4,4)+(n*seno**2)/long**2
           jacob(4,5)=jacob(4,5)-(n*seno*coseno)/long**2                          
           jacob(5,1)=jacob(5,1)+(n*seno*coseno)/long**2
           jacob(5,2)=jacob(5,2)-(n*coseno**2)/long**2     
           jacob(5,4)=jacob(5,4)-(n*seno*coseno)/long**2   
           jacob(5,5)=jacob(5,5)+(n*coseno**2)/long**2
c
c               "VER PAGINAS 25 Y 32 DEL ARTICULO: COMPUTER AND STRUCTURES"
c---------------------------------------------------------------------          
c
        else if(opcion(2).ne."small")then
                stop
     1                "error en la matriz opcion"
        end if
c
c******************************************************************
c        CALCULO DE LA CONTRIBUCION DINAMICA AL
c        JACOBIANO  GLOBAL
c******************************************************************
c
        if (opcion(1).eq."dynamic")then
c
                beta=param(2)
c
                call cal_masa(prop,max_prop,max_coo,max_nud,
     e                  max_desp,max_def,coord,
     s                  masa)
c     
                call cal_deriva_u(beta,dtiempo,max_for,
     s                        deriva_u)
c
                call pro_2mat(max_desp,max_desp,
     e                        max_desp,max_desp,
     e                        max_desp,max_desp,
     e                        max_desp,max_desp,max_desp,
     e                        masa,deriva_u,
     s                        jacob_iner)
c
                do i=1,nq
                        do j=1,nq
                           jacob(i,j)=jacob(i,j)+jacob_iner(i,j)
                        end do
                end do
c
        else if(opcion(1).ne."static")then
                stop
     1                        "error en la matriz opcion"
        end if
c
        return
        end
c
c
c******************************************************************
c
c        SUBRUTINA CAL_TANG
c        CALCULO DEL JACOBIANO LOCAL
c        
c        EL JACOBIANO SE CALCULA RESOLVIENDO EL SIGUIENTE
c        SISTEMA DE ECUACIONES MATRICIALES LINEALES;
c        
c        
c        |  [@R/@esf]   [@R/@INT_V]   [@R/@INT_A] | |  [@esf/@def]   |     | [@R/@def]  |
c        |                                        | |                |     |            |
c        |  [@T/@esf]   [@T/@INT_V]   [@R/@INT_A] | | [@INT_V/@def]  | = - |    [0]     |
c        |                                        | |                |     |            |
c        |  [@G/@esf]   [@G/@INT_V]   [@G/@INT_A] | | [@INT_A/@def]  |     |    [0]     |
c
c
c        DONDE
c
c       {R} representa a la ley de estado
c
c        Ti=  f(i)                                        si la variable interna i esta activa
c        Ti = varint_v(i) - varint_v_cero(i)                si la variable interna i no esta activa 
c
c        f(i) es la funcion de carga de la variable interna Vi
c
c        {G} define la Ley de Da~no
c
c        El simbolo @ representa la derivada parcial,  
c        ESF los esfuerzos        
c        INT_V las variables internas con funcin de fluencia
c        INT_A las variables internas sin funcin de fluencia
c
c        [@esf/@def]=tang
c        
c        ESCRITA POR MARÕA EUGENIA MARANTE y JULIO FLOREZ LOPEZ 
c        (Abril 2001)
c     Modificada por : Nayive Jaramillo S. (Enero 2002)
c******************************************************************
c
        subroutine cal_tang(max_prop,prop,max_jprop,jprop,max_opcion,
     e                     opcion,n_prop,max_var_v,max_var_a,
     e                     max_def,esf,varint_v_cero,varint_a_cero,
     e                     deftot,varint_v,varint_a,
     e                     esf_cero,deftot_cero,
     e                     param,max_param,max_act,activas,
     s                     tang)
c
c
       IMPLICIT NONE
c******************************************************************
c        DECLARACION DE VARIABLES MDC
c******************************************************************
c
c        prop= matriz de propiedades del elemento
c        param= parametros de la integracion numerica
c        param(1)=alfa
c        param(2)=beta
c        param(3)=gamma
c        param(4)=maximo numero de iteraciones
c        param(5)=maximo numero de pasos de integracin = PNEWDT
c        param(7:6+nfi+nv+na)=ceros numâricos de {R},{T},{G}
c        opcion= opciones del calculo (dinamico o estatico,
c                peque~os o grandes desplazamientos, etc)
c        opcion(1)= "static" o "dynamic"
c        opcion(2)= "small" o "large"
c        dtiempo= incremento de tiempo en el paso
c        varint_v_cero= variables internas con funcin de fluencia al principio del paso
c        varint_v = variables internas con funcin de fluencia al final del paso
c        varint_a_cero= variables internas sin funcin de fluencia al principio del paso
c        varint_a = variables internas sin funcin de fluencia al final del paso
c        nv= nîmero de variables internas con funcin de fluencia
c        na= nîmero de variables internas sin funcin de fluencia
c        nfi= dimensin de la matriz de deformaciones y de esfuerzos
c        n_prop= numero de propiedades del elemento
c        max_prop= tama~o de la matriz de propiedades
c        max_var_v= tama~o de la matriz de variables internas con funcin de fluencia
c        max_var_a= tama~o de la matriz de variables internas sin funcin de fluencia
c        max_opcion= tama~o de la matriz de opciones
c        max_param= tama~o de la matriz de parametros
c        max_def= tama~o de la matriz de deformaciones
c        esf=esfuerzos al final del paso
c        deftot=deformaciones totales al final del paso
c        activas= matriz que indica el estado de las variables internas
c
c******************************************************************
c
        integer max_prop,max_var_v,max_var_a,max_opcion,max_param
c
        integer max_jprop,max_def,max_act
c
        Real(kind=8):: prop(max_prop),param(max_param)
c
        character(len=8):: opcion(max_opcion)
c
        integer n_prop,jprop(max_jprop)
c
        Real(kind=8):: varint_v(max_var_v),varint_v_cero(max_var_v)
c
        Real(kind=8):: varint_a(max_var_a),varint_a_cero(max_var_a)
c
        Real(kind=8):: deftot(max_def),esf(max_def)
c        
        Real(kind=8):: esf_cero(max_def),deftot_cero(max_def)
c
        Real(kind=8):: tang(max_def,max_def)
c
        logical activas(max_act)
c
c******************************************************************
c                DECLARACIONES LOCALES
c******************************************************************
c
c        hiper= hipermatriz para la resolucion del sistema de ecuaciones
c                   locales por el metodo de Newton
c        max_hiper= tama~o maximo de hiper
c        max_esf= dimension maxima de deriv_r_def=max_def
c        deriv_r_def= matriz de derivadas de {R} con respecto a DEF
c        ind= matriz de almacenamiento de los terminos independientes
c
c******************************************************************
c
        integer max_hiper,max_esf,nfi,nv,na,nt,i,j,k,l
        parameter(max_hiper=11)
        parameter(max_esf=3)
c
        Real(kind=8):: hiper(max_hiper,max_hiper)
        Real(kind=8):: hiper_aux(max_hiper,max_hiper)
        Real(kind=8):: deriv_r_def(max_esf,max_esf)
        Real(kind=8):: ind(max_hiper)
        logical conv
c
c******************************************************************
c        CALCULO DE LA HIPERMATRIZ COMPUESTA POR LAS DERIVADAS DE
c        {R},{T} y {G} CON RESPECTO A esf,int_v e int_a
c******************************************************************
c
c
          call cal_hiper(max_prop,prop,max_opcion,opcion,
     e            max_jprop,jprop,n_prop,
     e            max_var_v,max_var_a,max_param,param,
     e            varint_v_cero,varint_a_cero,
     e            max_def,deftot_cero,deftot,esf_cero,
     e            varint_v,varint_a,esf,max_hiper,
     e            activas,max_act,
     s            hiper,conv)
c
c******************************************************************
c
          nfi=jprop(4)
          nv=jprop(5)
          na=jprop(6)
          nt=nfi+nv+na 
c
c*********************************************************************
c               LLAMADA A LA SUBRUTINA QUE CALCULA LA
c               DERIVADA DE R CON RESPECTO A DEF
c*********************************************************************
c
        call cal_deriv_r_def(max_prop,prop,max_opcion,opcion,
     e         max_jprop,jprop,n_prop,max_esf,
     e         max_var_v,max_var_a,max_param,param,varint_v_cero,
     e         varint_a_cero,max_def,
     e         deftot_cero,deftot,esf_cero,
     e         varint_v,varint_a,esf,max_hiper,
     s         deriv_r_def)
c
c         
        do j=1,nfi
c
           do k=1,nt
               do l=1,nt
                   hiper_aux(k,l)=hiper(k,l)
               end do
           end do
c
                do i =1,nfi
                        ind(i)=deriv_r_def(i,j)
                end do
c
                do i=nfi+1,nt
                        ind(i)=0.
                end do
c
c********************************************************************
c       RESOLUCION DEL SISTEMA DE ECUACIONES
c********************************************************************
c
          call rfssedp(hiper_aux,ind,nt)
c
c********************************************************************
c
                do i=1,nfi
                        tang(i,j)=-ind(i)
                end do
c
        end do
c
        return
        end
c
c
c
c******************************************************************
c
c        SUBRUTINA  CAL_DERIV_R_DEF
c        CALCULO DE LAS DERIVADAS DE {R} CON RESPECTO A {DEF}
c        DONDE {R}={f-fp}-[F(D+)]{M+}-[F(D-)]{M-}
c
c        ESCRITA POR J. FLOREZ LOPEZ y MARÕA EUGENIA MARANTE
c        (Mayo 2001)
c        Modificada por: Nayive Jaramillo S.(ENERO 2002)
c******************************************************************
c
        subroutine cal_deriv_r_def(max_prop,prop,max_opcion,opcion,
     e            max_jprop,jprop,n_prop,max_esf,
     e            max_var_v,max_var_a,max_param,param,varint_v_cero,
     e            varint_a_cero,max_def,
     e            deftot_cero,deftot,esf_cero,
     e            varint_v,varint_a,esf,max_hiper,
     s            deriv_r_def)
c
       IMPLICIT NONE
c******************************************************************
c                DECLARACION DE VARIABLES NUCLEO
c******************************************************************
c
c        prop= matriz de propiedades del elemento
c        param= parametros de la integracion numerica
c        param(1)=alfa
c        param(2)=beta
c        param(3)=gamma
c        param(4)=mﬂximo nînero de iteraciones
c        param(5)=mﬂximo nîmero de pasos de integracin =PNEWDT
c        param(7:6+nt)=ceros numericos para {R},{T},{X}
c        opcion= opciones del calculo (dinamico o estatico,
c                peque±os o grandes desplazamientos, etc)
c        opcion(1)= "static" o "dynamic"
c        opcion(2)= "small" o "large"
c        dtiempo= incremento de tiempo en el paso
c        varint_v_cero= variables internas con funcin de fluencia al principio del paso
c        varint_a_cero= variables internas sin funcin de fluencia al principio del paso
c        varint_v= variables internas con funcin de fluencia al final del paso
c        varint_a= variables internas sin funcin de fluencia al final del paso
c        nv= nîmero de variables internas con funcin de fluencia
c        na= nîmero de variables internas sin funcin de fluencia
c        ny= nîmero de fuerzas termodinﬂmicas
c        n_prop= nîmero de propiedades del elemento
c        nfi= dimensin de la matriz de deformaciones y de esfuerzos
c        max_prop= tama±o de la matriz de propiedades
c        max_var_v= tama±o de la matriz de var.int. con funcin de fluencia
c        max_var_a= tama±o de la matriz de var.int. sin funcin de fluencia 
c        max_opcion= tama±o de la matriz de opciones
c        max_param= tama±o de la matriz de parametros
c        max_def= tama±o de la matriz de deformaciones
c        esf=esfuerzos al final del paso
c        deftot=deformaciones totales al final del paso
c        esf_cero= esfuerzos al principio del paso
c        deftot_cero= deformaciones totales al principio del paso
c        fterm= Fuerzas termodinamicas asociadas a las variables 
c                    internas al final del paso
c        activas= matriz que indica el estado de todas las variables
c                  internas
c        max_act= tama±o de la matriz activas
c        rtg=({R},{T},{G})
c
c******************************************************************
c
        integer max_prop,max_var_v,max_var_a,max_opcion,max_param
c
        integer max_def,max_hiper,max_esf,max_jprop
c
        Real(kind=8):: prop(max_prop),param(max_param)
c
        character(len=8):: opcion(max_opcion)
c
        integer n_prop,jprop(max_jprop)
c
        Real(kind=8):: varint_v(max_var_v),varint_a(max_var_a),
     2           varint_v_cero(max_var_v),varint_a_cero(max_var_a)
c
        Real(kind=8):: deftot_cero(max_def),esf_cero(max_def),
     1                deftot(max_def),esf(max_def)
c
        Real(kind=8):: deriv_r_def(max_esf,max_esf)
c
c******************************************************************
c
        deriv_r_def(1,1)=1.
        deriv_r_def(1,2)=0.
        deriv_r_def(1,3)=0.
        deriv_r_def(2,1)=0.
        deriv_r_def(2,2)=1.
        deriv_r_def(2,3)=0.
        deriv_r_def(3,1)=0.
        deriv_r_def(3,2)=0.
        deriv_r_def(3,3)=1.
c
        return
        end
c
c******************************************************************
c******************************************************************
c
c        SUBRUTINA cal_masa
c        CALCULO DE LA MATRIZ DE MASA
c        MODELO DE MASA CONCENTRADA
c
c        ESCRITA POR MARIA EUGENIA MARANTE y JULIO FLOREZ LOPEZ
c     (mayo 2001)
c     MODIFICADA POR: NAYIVE JARAMILLO (ENERO 2002)
c******************************************************************
c
        subroutine cal_masa(prop,max_prop,max_coo,max_nud,
     e                  max_desp,max_def,coord,
     s                  masa)
c     
c
       IMPLICIT NONE
c******************************************************************
c                DECLARACION DE VARIABLES 
c******************************************************************
c
c        prop= matriz de propiedades del elemento
c        max_prop= tama~o de la matriz de propiedades
c        max_desp=maximo numero de despalazamientos
c        max_def=tama~o de matriz de deformaciones
c        masa=matriz de masa consistente global
c
c******************************************************************
c
        integer max_desp,max_prop,max_def,max_coo,max_nud
c
        Real(kind=8):: masa(max_desp,max_desp)
c
        Real(kind=8):: prop(max_prop)
        
        Real(kind=8):: coord(max_coo,max_nud)
c
c******************************************************************
c                DECLARACIONES LOCALES
c******************************************************************
c        
c        masa_total=masa total del elemento
c        l=longitud del elemento
c        s=seno , c=coseno
c
c******************************************************************
c
        Real(kind=8):: masa_total,l,s,c
c
c******************************************************************
c              
         l=prop(1)
         s=prop(2)
         c= prop (3)
         masa_total=prop(22)
c
c******************************************************************
c
        masa(1,1)=(35.*c*c+39.*s*s)*masa_total/105.
        masa(1,2)=-16.*c*s*masa_total/420.
        masa(1,3)=22.*l*s*masa_total/420.
        masa(1,4)=(35.*c*c+27.*s*s)*masa_total/210.
        masa(1,5)=16.*c*s*masa_total/420.
        masa(1,6)=-13.*l*s*masa_total/420.
        masa(2,1)=masa(1,2)
        masa(2,2)=(39.*c*c+35.*s*s)*masa_total/105.
        masa(2,3)=-22.*c*l*masa_total/420.
        masa(2,4)=16.*c*s*masa_total/420.
        masa(2,5)=(27.*c*c+35.*s*s)*masa_total/210.
        masa(2,6)=13.*c*l*masa_total/420.
        masa(3,1)=masa(1,3)
        masa(3,2)=masa(2,3)
        masa(3,3)=l*l*masa_total/105.
        masa(3,4)=13.*l*s*masa_total/420.
        masa(3,5)=-13.*c*l*masa_total/420.
        masa(3,6)=-l*l*masa_total/140.
        masa(4,1)=masa(1,4)
        masa(4,2)=masa(2,4)
        masa(4,3)=masa(3,4)
        masa(4,4)=masa(1,1)
        masa(4,5)=masa(1,2)
        masa(4,6)=-masa(1,3)
        masa(5,1)=masa(1,5)
        masa(5,2)=masa(2,5)
        masa(5,3)=masa(3,5)
        masa(5,4)=masa(4,5)
        masa(5,5)=masa(2,2)
        masa(5,6)=-masa(2,3)
        masa(6,1)=masa(1,6)
        masa(6,2)=masa(2,6)
        masa(6,3)=masa(3,6)
        masa(6,4)=masa(4,6)
        masa(6,5)=masa(5,6)
        masa(6,6)=masa(3,3)
c
        return
        end
c
c
c******************************************************************
c
c                SUBRUTINA cal_deriva_u
c                CALCULO DE LA MATRIZ DE DERIVADAS DE LA ACELERACION
c                CON RESPECTO A LOS DESPLAZAMIENTOS
c
c                ENTRADA: LA PRIMERA LINEA DE LA SUBRUTINA
c                SALIDA LA SEGUNDA LINEA
c
c                VERSION 3.0 (06/10/93)
c
c******************************************************************
c
        subroutine cal_deriva_u(beta,dtiempo,max_for,
     s                        deriva_u)

c
       IMPLICIT NONE
c******************************************************************
c        
c        deriva_u= derivadas de las aceleraciones con respecto a los
c                        desplazamientos
c        beta= parametro de la integracion numerica en el caso dinamico
c        dtime= incremento de tiempo
c
c******************************************************************
c
        integer max_for
c
        Real(kind=8):: beta,dtiempo
        Real(kind=8):: deriva_u(max_for,max_for)
        integer i,j
c *****************************************************************
        do i=1,max_for
            do j=1,max_for
                 deriva_u(i,j)=0.0
             end do
        end do                
c********************************************************************
        if(dtiempo.eq.0.)then
            do i=1,max_for
                deriva_u(i,i)=0.
            end do
        else        
            do i=1,max_for
                deriva_u(i,i)=1./(beta*dtiempo*dtiempo)
            end do
        end if
c
        return
        end
c
c
c
c
c******************************************************************************
c
c       CALCULO DE MOMENTO DE AGRIETAMIENTO (MCR),
c       MOMENTO PLASTICO (MP),
c       MOMENTO ULTIMO (MU) Y 
c       CURVATURA ULTIMA (XU)
c
c       ESCRITO POR: RICARDO PIC”N, SCARLET Y JULIO FLOREZ
c       INICIADA EL: 28 DE MAYO DE 2021
c
c******************************************************************************
c
       subroutine cal_generador_1(
     e             max_prop,prop,na,ei,l,esf,
     s             gcr_pos_i,q_pos_i,my_pos_i,c_pos_i,dp_pos_i,du_pos_i,
     s             gcr_neg_i,q_neg_i,my_neg_i,c_neg_i,dp_neg_i,du_neg_i,
     s             gcr_pos_j,q_pos_j,my_pos_j,c_pos_j,dp_pos_j,du_pos_j,
     s             gcr_neg_j,q_neg_j,my_neg_j,c_neg_j,dp_neg_j,du_neg_j)
c
       IMPLICIT NONE
c******************************************************************************
c       DECLARACIONES
c
c
c       ACEROS PARA EL EXTREMO i
c
c       as_posit_i=area de acero positiva (area de acero en traccion
c               cuando el momento es positivo en el extremo i del
c               elemento)
c       as_negat_i=area de acero negativa en el extremo i
c       rec_posit=recubrimiento positivo en el extremo i
c       rec_negat=recubrimiento negativo en el extremo i
c
c       ACEROS PARA EL EXTREMO J
c
c       as_posit_j=area de acero positiva (area de acero en traccion
c               cuando el momento es positivo en el extremo j del
c               elemento)
c       as_negat_j=area de acero negativa en el extremo j
c       rec_posit_j=recubrimiento positivo en el extremo j
c       rec_negat_j=recubrimiento negativo en el extremo j
c
c       na: fuerza axial
c       
c*******unidades*******unidades*****unidades*****unidades*****unidades*****
c
c		Nota: Las unidades de entrada son kilogramo fuerza y centimetros
c
c			  Las unidades de salida son toneladas fuerza y centimetros
c*******unidades*******unidades*****unidades*****unidades*****unidades*****
c***************************************************************************
	integer max_prop,i,ng,max_props
	parameter(max_props=50)
c
	real*8 props(max_props),esf(3)
	real*8 prop(max_prop),as_posit_i,as_negat_i,rec_posit_i
	real*8 rec_negat_i,as_posit_j,as_negat_j,rec_posit_j,rec_negat_j
	real*8 as1_i,d1_i,as2_i,d2_i,as3_i,d3_i,as4_i,d4_i,as5_i,d5_i
	real*8 as1_in,d1_in,as2_in,d2_in,as3_in,d3_in,as4_in,d4_in,as5_in
	real*8 d5_in
	real*8 as1_j,d1_j,as2_j,d2_j,as3_j,d3_j,as4_j,d4_j,as5_j,d5_j
	real*8 as1_jn,d1_jn,as2_jn,d2_jn,as3_jn,d3_jn,as4_jn,d4_jn,as5_jn
	real*8 d5_jn
	real*8 ei,ea,masa,ei_mod,ci,cj,na
c
	real*8 b,h,fccmin
	real*8 fy,fsu,ey,esh,esm,fyh,mesy,mesh
	real*8 hc,zm,k,ps,sh,rh,rv,bc
	real*8 l,eccu_u
	real*8 fc,eo,euc,eum,mec,eccu
	character*14 op_diseno,op_deform
	character*14 extremo_libre
	real*8 E1,a1,a2
	real*8 ao,de,rx,ry,le
      Real*8 gcr_pos_i,q_pos_i,my_pos_i,c_pos_i
	Real(kind=8):: dp_pos_i,du_pos_i
      Real(kind=8):: gcr_neg_i,q_neg_i,my_neg_i
      Real(kind=8):: c_neg_i,dp_neg_i,du_neg_i
	Real(kind=8):: gcr_pos_j,q_pos_j,my_pos_j
      Real(kind=8):: c_pos_j,dp_pos_j,du_pos_j
	Real(kind=8):: gcr_neg_j,q_neg_j,my_neg_j
	Real(kind=8):: c_neg_j,dp_neg_j,du_neg_j
c
c******************************************************************************
c       LECTURA DE DATOS
c******************************************************************************
c
      as_posit_i=prop(17)
	as_negat_i=prop(19)
	rec_posit_i=prop(18)
      rec_negat_i=prop(20)
	masa=prop(38)
      as_posit_j=prop(41)
	as_negat_j=prop(43)
	rec_posit_j=prop(42)
	rec_negat_j=prop(44)
        as1_i=prop(21)
	d1_i=prop(22)
	as2_i=prop(23)
	d2_i=prop(24)
	as3_i=prop(25)
	d3_i=prop(26)
	as4_i=prop(27)
	d4_i=prop(28)
	as5_i=prop(29)
      d5_i=prop(30)
      as1_j=prop(45)
	d1_j=prop(46)
	as2_j=prop(47)
	d2_j=prop(48)
	as3_j=prop(49)
	d3_j=prop(50)
	as4_j=prop(51)
	d4_j=prop(52)
	as5_j=prop(53)
      d5_j=prop(54)
c
c
c
	b=prop(4)
	h=prop(5)
      fc=prop(7)
	eo=prop(8)
	euc=prop(9)
      fy=prop(11)
	fsu=prop(12)
	ey=prop(13)
	esh=prop(14)
	esm=prop(15)
	fyh=prop(16)
	de=prop(31)
      rx=prop(32)
	ry=prop(33)
	le=prop(34)
	rh=prop(35)
	rv=prop(36)
	sh=prop(37)
      E1=prop(39)
	eccu_u=prop(40)
c
c	
c	
	props(28)=as1_j+as2_j+as3_j+as4_j+as5_j
	props(29)=as1_i+as2_i+as3_i+as4_i+as5_i
	props(1)=b
	props(2)=h
	props(30)=l
	props(31)=rx
	props(32)=ry
	props(33)=de
	props(34)=le
	props(35)=eccu_u
c
c       OpciÛn de DiseÒo (Confinado=1  Û   No Confinado=0)  ==>  prop(6)
c
	if ((prop(6).eq.0.).or.(prop(6).eq.1.)) then
	     if (prop(6).eq.1.) then
		 props(3)=1.
c		C·lculo del ancho del nucleo confinado
		 bc=b-2.*rh
		 hc=h-2.*rv
c
     		 ao=(de*de*3.141592654)/4.
		 if ((rx.eq.0.).and.(ry.eq.0.).and.(le.gt.0.)) then
			ps=ao*le/(sh*bc*hc)
		 else if ((rx.ge.2.).and.(ry.ge.2.).and.(le.eq.0.)) then
			ps=ao*(rx*bc+ry*hc)/(sh*bc*hc)
		 end if
c	  CondiciÛn BAJA significa que las cargas son aplicadas LENTAMENTE 
c	  (Est·ticas)
c
c       OpciÛn de Velocidad de deformaciÛn (Baja=1  Û   Alta=0)  ==>  prop(10)
c
		if(prop(10).eq.1.) then
			k=(1.+(ps*fyh)/fc)
			zm=0.5/((30.6+0.29*fc)/(145.*fc-10200.)+
     1			   (0.75*ps*DSQRT(hc/sh))-k*0.002)
			eccu=(0.8+zm*k*0.002)/zm
			eum=0.004+(0.9*ps*fyh)/3060.
c       C·lculo de fccmin porque se hace siempre con eccu calculado
			a1=0.2*k*fc
			a2=k*fc*(1.-0.8*((eum-0.002*k)/(eccu-0.002*k)))
			fccmin=MAX(a1,a2)
c
			if (eccu.ge.eum) then
				eccu=eum
			end if
c
c	  CondiciÛn ALTAS significa que las cargas son aplicadas RAPIDAMENTE 
c	  (Din·micas)
		else  ! deformaciones altas
c
			k=1.25*(1.+(ps*fyh)/fc)
			zm=1.25*0.5/((30.6+0.29*fc)/(145.*fc-10200.)+
     1                     (0.75*ps*DSQRT(hc/sh))-k*0.002)
			eccu=(0.8+zm*k*0.002)/zm
			eum=0.004+(0.9*ps*fyh)/3060.
c
c             C·lculo de fccmin porque se hace siempre con eccu calculado
			a1=0.2*k*fc
			a2=k*fc*(1.-0.8*((eum-0.002*k)/(eccu-0.002*k)))
			fccmin=MAX(a1,a2)
c
			if (eccu.ge.eum) then
				eccu=eum
			end if
		end if
c***********************************************************************************
c***********************************************************************************
c
c			Deformacion ultima en el concreto confinado Para el analisis
c
c***********************************************************************************
c***********************************************************************************
c
		if ((eccu.le.eccu_u).and.(eccu_u.le.eum)) then
			eum=eccu_u
		end if
c
		if ((eccu_u.le.eccu).and.(eccu.le.eum).and.
     1               (eccu_u.gt.0.)) then
			eum=eccu_u
			eccu=eccu_u
		end if
c
		if (eccu_u.gt.eum) then
		   write(*,*)"****	****		AVISO	****	*****"
		   write(*,*)"Fisicamente la deformacion del concreto confinado"
		   write(*,*)"en este analisis no puede ser mayor a :",eum
		   write(*,*)"El programa tomara como deformacion "
		   write(*,*)"maxima en el concreto confinado el valor anterior"
		end if
c
	     props(6)=euc
		 props(12)=fyh
		 props(13)=ps
		 props(14)=rh
		 props(15)=rv
		 props(16)=sh
	     props(26)=eccu
		 props(27)=fccmin
c
	     end if 
	     if (prop(6).eq.0.) then
				props(3)=0.
				props(6)=euc
				props(12)=fyh
				props(14)=rh
				props(15)=rv
				props(16)=sh
           end if
c
	else
	    stop "error en opcion de diseno"
c	
	end if
c
	props(4)=fc
	props(5)=eo
	props(7)=fy
	props(8)=fsu
	props(9)=ey
	props(10)=esh
	props(11)=esm
c
	mesy=fy/ey
	props(17)=mesy
	mesh=2.*(fsu-fy)/(esm-esh)
	props(18)=mesh
	mec=E1
	props(19)=mec
	props(20)=k
	props(21)=hc
	props(22)=eum
	props(23)=zm
	props(24)=0.
	props(25)=0.
c
c	Distancias a cada nivel de acero en las condiciones negativas
c	CondiciÛn i negativa
c
	d1_in=h-d5_i
	d2_in=h-d4_i
	d3_in=h-d3_i
	d4_in=h-d2_i
	d5_in=h-d1_i
	as1_in=as5_i
	as2_in=as4_i
	as3_in=as3_i
	as4_in=as2_i
	as5_in=as1_i
c	CondiciÛn J negativa
	d1_jn=h-d5_j
	d2_jn=h-d4_j
	d3_jn=h-d3_j
	d4_jn=h-d2_j
	d5_jn=h-d1_j
	as1_jn=as5_j
	as2_jn=as4_j
	as3_jn=as3_j
	as4_jn=as4_j
	as5_jn=as1_j
c
c
c******************************************************************************
c       DETERMINACION DE LOS DIAGRAMAS DE INTERACCION PARA
c       MOMENTO POSITIVO EN EL EXTREMO i
c******************************************************************************
c
      call diagramas_f(max_props,props,as_posit_i,rec_posit_i,esf,na,l,
     e		   as1_i,d1_i,as2_i,d2_i,as3_i,d3_i,as4_i,d4_i,1,
     e             as5_i,d5_i,as_negat_i,rec_negat_i,
     s             gcr_pos_i,q_pos_i,my_pos_i,c_pos_i,
     s             dp_pos_i,du_pos_i)
c
c******************************************************************************
c
c       DETERMINACION DE LOS DIAGRAMAS DE INTERACCION PARA
c       MOMENTO NEGATIVO EN EL EXTREMO i
c******************************************************************************
c
      call diagramas_f(max_props,props,as_negat_i,rec_negat_i,esf,na,l,
     e	      as1_in,d1_in,as2_in,d2_in,as3_in,d3_in,as4_in,d4_in,3,
     e        as5_in,d5_in,as_posit_i,rec_posit_i,
     s        gcr_neg_i,q_neg_i,my_neg_i,c_neg_i,
     s        dp_neg_i,du_neg_i)
c
c
c******************************************************************************
c       DETERMINACION DE LOS DIAGRAMAS DE INTERACCION PARA
c       MOMENTO POSITIVO EN EL EXTREMO J
c******************************************************************************
c
      call diagramas_f(max_props,props,as_posit_j,rec_posit_j,esf,na,l,
     e 	             as1_j,d1_j,as2_j,d2_j,as3_j,d3_j,as4_j,d4_j,2,
     e               as5_j,d5_j,as_negat_j,rec_negat_j,
     s               gcr_pos_j,q_pos_j,my_pos_j,c_pos_j,
     s               dp_pos_j,du_pos_j)
c
c******************************************************************************
c
c       DETERMINACION DE LOS DIAGRAMAS DE INTERACCION PARA
c       MOMENTO NEGATIVO EN EL EXTREMO J
c******************************************************************************
c
      call diagramas_f(max_props,props,as_negat_j,rec_negat_j,esf,na,l,
     e	      as1_jn,d1_jn,as2_jn,d2_jn,as3_jn,d3_jn,as4_jn,d4_jn,4,
     e        as5_jn,d5_jn,as_posit_j,rec_posit_j,
     s        gcr_neg_j,q_neg_j,my_neg_j,c_neg_j,
     s        dp_neg_j,du_neg_j)
c
c
c******************************************************************************
c
c
	RETURN
	END
c	
c******************************************************************************
c
c******************************************************************************
c       DIAGRAMA: CALCULO DE LOS DIAGRAMAS DE INTERACCION
c       ESCRITA POR MARIA ELENA PERDOMO Y JULIO FLOREZ LOPEZ
c       COMENZADA: 21 DE JULIO DE 1997
c******************************************************************************
c
	subroutine diagramas_f(max_props,props,as_trac,rec_trac,esf,na,l,
     e				as1,d1,as2,d2,as3,d3,as4,d4,rotula,
     e				as5,d5,as_comp,rec_comp,
     s                gcr,q,m,c,
     s                dp,du)
c
       IMPLICIT NONE
c******************************************************************************
c       DECLARACIONES GLOBALES
c******************************************************************************
c
	integer max_props,rotula
	real*8 props(max_props),as_trac,as_comp,rec_trac,rec_comp,esf(3)
	real*8 as1,d1,as2,d2,as3,d3,as4,d4,as5,d5
	real*8 tabla_mcr(4,2)
	real*8 tabla_mp(4,3)
	real*8 tabla_mu(5,3)
	real*8 tabla_x(4,2)
	real*8 ei_mod,mcr,mp,mu,x,na,dp,du,l,gcr,q,m,c
c
c******************************************************************************
c       DECLARACIONES LOCALES
c******************************************************************************
c
	real*8 b,h
	real*8 hc
	real*8 fccmin,eccu
	real*8 fy,fsu,ey,esh,esm,fyh,mesy,mesh
     	real*8 zm,k,ps,sh,rh,rv
	real*8 fc,eo,euc,eum,mec
	real*8 fr
      logical confinado
	real*8 area,inerc
	real*8 inercia_sg
	real*8 d,y,yt
	real*8 ea,ei,gav,av,g
	character*14 op_diseno
	real*8 n
	integer w
	parameter (w=5)
	real*8 di1(w),as_interm(w)
	real*8 zi,zj,z
	real*8 mi,mj,lp
c
c******************************************************************************
c
	b=props(1)
	h=props(2)
	if(props(3).eq.1)then
		confinado=.true.
		op_diseno="confinado"
	else ! props(3)=0
		confinado=.false.
		op_diseno="noconfinado"
	end if
c
     	eccu=props(26)
	fccmin=props(27)
	euc=props(6)
	fyh=props(12)
	ps=props(13)
	rh=props(14)
	rv=props(15)
	sh=props(16)
	fc=props(4)
	eo=props(5)
	fy=props(7)
	fsu=props(8)
	ey=props(9)
	esh=props(10)
	esm=props(11)
c
	mesy=props(17)
	mesh=props(18)
	mec=props(19)
	k=props(20)
	hc=props(21)
	eum=props(22)
      zm=props(23)
c       Calculo del area gruesa de la seccion
	area=b*h
c
c        Calculo de la distancia del eje centroidal a la fibra
c        m·s traccionada
c
       	n=mesy/mec
	 yt=(area*h/2.+(n-1.)*(as1*(h-d1)+as2*(h-d2)+as3*(h-d3)+
     1       as4*(h-d4)+as5*(h-d5)+
     1       (as_trac*rec_trac+as_comp*(h-rec_comp))))/
     1       (area+(n-1.)*(as_trac+as_comp+as1+as2+as3+as4+as5))
c
c             Calculo de la inercia de la seccion con el acero
c
	inerc=b*(h*h*h)/12.+area*(h/2.-yt)*(h/2.-yt)+(n-1.)*
     1     (as_trac*(yt-rec_trac)*(yt-rec_trac)+as_comp*
     1     (h-yt-rec_comp)*(h-yt-rec_comp))+(n-1.)*(as1*(yt-d1)*(yt-d1)+
     1      as2*(yt-d2)*(yt-d2)+as3*(yt-d3)*(yt-d3)+as4*(yt-d4)*(yt-d4)+
     1      as5*(yt-d5)*(yt-d5))
c
c			Calculo de la inercia de la secciÛn gruesa
c
	inercia_sg=b*h*h*h/12.
c
c       Calculo del modulo de elasticidad*inercia seccion transformada (se 
c	  divide entre 1000. para llevarlo a ton*cm^2)
c
	ei=mec*inerc/1000.
c
c       Calculo del modulo de elasticidad*area (se divide entre 1000. para 
c       llevarlo a ton*cm^2)
c
	ea=mec*area/1000.
c
	ei_mod=mec*inercia_sg/1000.
c
c*********************************************************************************
c       Calculo del modulo de corte (G) *area de corte efectiva (Av) 
c       (se divide entre 1000. para llevarlo a ton)
c******************************************************************************
c
      av=0.8*area
	g=mec/2./1.2
	gav=g*av/1000.
c
c
c******************************************************************************
c******************************************************************************
c       Calculo del modulo de rotura del concreto
	fr=2.*DSQRT(fc)
c	print*,fr
c
c******************************************************************************
c     Calculo de la altura util
	d=h-rec_trac
c
	y=h-yt
c	Distancias dedes la fibra m·s comprimida a c/nivel de acero
	di1(1)=d1
	di1(2)=d2
	di1(3)=d3
	di1(4)=d4
	di1(5)=d5
c     Vector de las areas de los aceros intermedios
	as_interm(1)=as1
	as_interm(2)=as2
	as_interm(3)=as3
	as_interm(4)=as4
	as_interm(5)=as5
c
      mi=esf(1)
	mj=esf(2)
c
c***************************************************************************
c
	call cal_mcr(fr,inerc,yt,h,
     e             mesy,as_trac,as_comp,as1,as2,as3,as4,as5,b,
     e             fc,mec,
     s             tabla_mcr)
c
      call interpolacion(na,4,tabla_mcr,4,
     s                   mcr)
c
	call cal_mp(op_diseno,
     e             h,fy,fsu,ey,esh,esm,
     e             mesy,mesh,as_trac,as_comp,rec_trac,rec_comp,w,
     e			 as_interm,di1,d,b,
     e             rv,rh,ps,sh,fyh,k,zm,fc,eo,euc,eum,eccu,
     e             hc,fccmin,
     s             tabla_mp,tabla_mcr)
c
      call interpolacion(na,4,tabla_mp,4,
     s	               mp)
c
c
	call cal_mu(op_diseno,tabla_mp,
     e        h,fy,fsu,ey,esh,esm,
     e        mesy,mesh,as_trac,as_comp,w,as_interm,di1,d,rec_trac,b,
     e		rec_comp,rv,rh,ps,sh,fyh,k,zm,fc,eo,euc,eum,eccu,
     e        hc,fccmin,
     s        tabla_mu)
c
c
      call interpolacion(na,5,tabla_mu,5,
     s                  mu)
c
	call cal_x(tabla_mp,tabla_mu,
     s                  tabla_x)
c
c
      call interpolacion(na,4,tabla_x,4,
     s                   x)
c	
c
      call cal_coefic(mcr,mp,mu,x,ei_mod,l,
     s                gcr,q,m,c,dp,du)
c
c******************************************************************
c       Calculo de la rotaciones plasticas ultimas
c******************************************************************
c
        mi=esf(1)
        mj=esf(2)
        if((mi.gt.0.and.mj.gt.0).or.(mi.lt.0.and.mj.lt.0))then
                zi=l*mi/(mi+mj)
                zj=l*mj/(mi+mj)
        else if(mi.eq.0.and.mj.eq.0)then
                zi=0.
                zj=0.
        else if(dabs(mi).gt.dabs(mj))then
                zi=l
                zj=0.
        else if(dabs(mi).lt.dabs(mj))then
                zi=0.
                zj=l
        else if((dabs(mi).eq.dabs(mj)).and.(mi.ne.mj))then
                zi=l
                zj=l
        else
                stop "error en el calculo de z"
        end if
c
        if ((rotula.eq.1).or.(rotula.eq.3)) then
	       z=zi
	  else
	       z=zj
	  end if
c
        lp=0.5*d+0.05*z
c
c*****************************************************************
c
        c=c/lp
c
c
c
c
	return
	end
c
c******************************************************************************
c       C¡LCULO DEL MOMENTO DE AGRIETAMIENTO
c
c       CAL_MCR
c       DIAGRAMA DE INTERACCION DEL DIAGRAMA DEL MOMENTO DE
c       AGRIETAMIENTO
c       ESCRITA POR MARIA ELENA PERDOMO Y JULIO FLOREZ LOPEZ
c       MODIFICADA POR RICARDO PIC”N, SCARLET Y JULIO FL”REZ 2021
c       INICIADA EL 21 DE JULIO DE 1997
c
c       modificada por: Nayive Jaramillo, 2002
c******************************************************************************
c
	subroutine cal_mcr(fr,inerc,y,h,
     e                   mesy,ast,asc,a1,a2,a3,a4,a5,b,
     e                   fc,mec,
     s                   tabla_mcr)
c
       IMPLICIT NONE
c******************************************************************************
c               DECLARACIONES LOCALES
c******************************************************************************
c
c       fr=modulo de rotura del concreto
c       inerc=inercia de la seccion
c       y=distancia del eje centroidal a la fibra mas traccionada
c       h=altura de la seccion
c       mesy=modulo de elasticidad del acero de refuerzo
c       mec=modulo de elasticidad del concreto
c       fc=resistencia maxima a compresion del concreto
c       mcr=momento de agrietamiento
c       area= area de la seccion transformada
c       b=ancho de la secci'on
c       ast=area de acero a traccion
c       asc=area de acero a compresion
c	  a1,a2,a3,a4,a5=areas de acero de las capas intermedias
c       na = Nivel de carga axial
c       c = Ìndice de corrosiÛn
c
c******************************************************************************
c
	real*8 fr,inerc,y,h,ast,asc,a1,a2,a3,a4,a5
	real*8 mesy,na,c
	real*8 fc,mec
	real*8 area !area transformada
	real*8 b
	real*8 tabla_mcr(4,2)
c
c****************************************************************************** 
c    LOS CUATRO PUNTOS DE MOMENTO AGRIETAMIENTO
c       Calaculo de la secciÛn transformada
c
	area=b*h+(mesy/mec-1.)*(ast+asc+a1+a2+a3+a4+a5)
c
c
	tabla_mcr(1,1)=-(-fr*area/1000.)
	tabla_mcr(1,2)=0.
	tabla_mcr(2,1)=0.
	tabla_mcr(2,2)=fr*inerc/(y*1000.)
	tabla_mcr(3,1)=-((y*(0.45*fc+fr)-h*fr)*area/(h*1000.))
	tabla_mcr(3,2)=(0.45*fc+fr)*inerc/(h*1000.)
	tabla_mcr(4,1)=-(0.45*fc*area/1000.)
	tabla_mcr(4,2)=0.
c
c
	return
	end
c
c****************************************************************************** 
c             C¡LCULO DEL MOMENTO PL¡STICO
c
c       CAL_MP
c       DIAGRAMA DE INTERACCION DEL DIAGRAMA DEL MOMENTO DE
c       PLASTICO
c       ESCRITA POR MARIA ELENA PERDOMO Y JULIO FLOREZ LOPEZ
c       MODIFICADA POR RICARDO PIC”N, SCARLET Y JULIO FL”REZ 2021
c       INICIADA EL 21 DE JULIO DE 1997
c
c	  modificada por: Nayive Jaramillo S., 2002
*******************************************************************************
c
	subroutine cal_mp(op_diseno,
     e                   h,fy,fsu,ey,esh,esm,
     e                   mesy,mesh,ast,asc,dt,dc,w,
     e			       as_interm,di1,d,b,
     e                   rv,rh,ps,sh,fyh,k,zm,fc,eo,euc,eum,eccu,
     e                   hc,fccmin,
     s                   tabla_mp,tabla_mcr)
c
       IMPLICIT NONE
c******************************************************************************
c               DECLARACIONES
c******************************************************************************
c
c       op_diseno=opciones de diseno para concreto confinado no
c                 confinado
c       fy=esfuerzo de cedencia del acero de refuerzo
c       fsu=esfuerzo maximo del acero de refuerzo , generalmente es 1.3Fy, y 
c	  puede lelgar a ser 1.6fy
c       ey=deformacion cedente del acero de refuerzo
c       esh=deformacion al final de la cedencia del acero
c       esm=deformacion maxima del acero de refuerzo
c       mesy=modulo de elasticidad del acero en la zona elastica
c       mesh=modulo de elasticidad del acero en la zona de endu-
c            recimiento.
c       ast=area de acero a traccion
c       asc=area de acero a compresion
c       d=altura util
c       dc=recubrimiento del acero a compresion
c       dt=recubrimiento del acero a traccion
c	  distancias desde la fibra m·s comprimida a c/capa intermedia de acero
c	  di1(1),di1(2),di1(3),di1(4),di1(5)
c       y=distancia del eje centroidal a la fibra mas traccio-
c         nada
c       h=altura de la seccion
c       hc=ancho del nucleo confinado
c       rv=recubrimiento vertical de los estribos
c       k=factor de incremento de resistencia del concreto con-
c         finado.
c       zm=parametro que define la pendiente de la rama descen-
c          diente
c       fc=resistencia maxima a compresion del concreto
c       eo=deformacion del concreto para la resistencia maxima
c       euc=deformacion ultima del concreto no comfinado o defor-
c           macion al final de la rama descendente del concreto
c           confinado
c       eum=deformacion ultima del concreto confinado
c       mp=momento plastico
c       x=profundidad del eje neutro
c       ec=deformacion del concreto en la fibra mas comprimida
c       eps=deformacion en el acero de refuerzo de la fibra  m·s comprimida
c       es=deformacion de la fibra mas traccionada
c       fs=esfuerzo del acero en la fibra mas traccionada
c       fps=esfuerzo del acero en la fibra m·s comprimida
c       cs=resultante de la fuerza de compresion del acero de
c          refuerzo
c       cc=resultante de las fuerzas de compresion del concreto
c       mc=momento resistente debido a las fuerzas de compresion
c          del concreto
c       trac=resultante de las fuerzas de traccion
c       comp=resultante de las fuerzas de compresion
c       fiy=rotacion de la seccion correspondiente a la primera
c           fluencia del acero.
c	  d11(i)=distancias al eje neutro, en cada capa de acero intermedia
c     fs_interm(i)= esfuerzo a nivel de los aceros intermedios
c     Fz_s_interm(i)= fuerza a nivel de los aceros intermedios
c     c_s= fuerza a compresiÛn a niveles de los aceros intermedios
c	cs1= fuerza a compresiÛn en el acero m·s comprimido
c	t= fuerza a tracciÛn a nivel de los aceros
c	m_s=momento resistente debido a las fuerzas de los aceros intermedios
c     aux1, aux2, correcciÛn de la fuerza de compresiÛn a nivel de los 
c	aceros
c     fs_interm(i)= esfuerzo a nivel de los aceros intermedios
c     f_comp(i)= esfuerzo a compresiÛn a nivel de los aceros intermedios
c     Fz_s_interm(i)= fuerza a nivel de los aceros intermedios
c     c_s= fuerza a compresiÛn a niveles de los aceros intermedios
c	cs1= fuerza a compresiÛn en el acero m·s comprimido
c	t= fuerza a tracciÛn a nivel de los aceros
c	m_s=momento resistente debido a las fuerzas de los aceros intermedios
c     aux1, aux2, correcciÛn de la fuerza de compresiÛn a nivel de los 
c	aceros
c	es_interm(i) deformaciÛn a nivel de los aceros intermedios
c******************************************************************************
c
	integer w
	character*14 op_diseno
	real*8 as_interm(w),di1(w)
	real*8 h,d,b
	real*8 ast,asc,dt,dc
	real*8 hc,ps,rv,rh,sh,k,zm
	real*8 fyh,c1,c2
	real*8 fy,fsu,ey,esh,esm
	real*8 fs,fps,eps
	real*8 mesy,mesh,fccmin
	real*8 fc,eo,euc,eum,ec,eccu
	real*8 tabla_mp(4,3)
	real*8 tabla_mcr(4,3)
c************************ Locales *********************************************
	real*8 x,cc,mc,cs,comp,trac,trac1
	real*8 f
	real*8 xini,xfin,Pplastico
	logical convergencia
	real*8 mp,pp
	integer i_ite,i
      real*8 as_t_inter
	real*8 d11(5),es_interm(5),fs_interm(5),f_compre(5)
	real*8 Fz_s_interm(5),c_s,m_s
	real*8 aux1,aux2,cs1,t,p
	real*8 Maux,Paux
c******************************************************************************
c       PRIMER VALOR DE LA TABLA
c       FLUENCIA DEL ACERO EN TRACCION Y SU CURVATURA CEDENTE                
c       PUNTO 1 (ABAJO) MOMENTO PLASTICO
c******************************************************************************
	as_t_inter=0.
	do i=1,5
		as_t_inter=as_interm(i)+as_t_inter
	end do
c
      p=-((ast+asc+as_t_inter)*fy)/1000.
	if (dabs(p).le.dabs(tabla_mcr(1,1))) then
	  tabla_mp(1,1)=-(1.05*tabla_mcr(1,1))
	else
        tabla_mp(1,1)=-p
      endif 
	tabla_mp(1,2)=0.
	tabla_mp(1,3)=0.
c******************************************************************************
c
c******************************************************************************
c       SEGUNDO PUNTO DE LA TABLA
c       CALCULO DEL MOMENTO DE FLUENCIA PARA CARGA AXIAL                 
c	  PUNTO 2 ( AXIAL = 0 ) MOMENTO PLASTICO
c       CERO Y SU CURVATURA CEDENTE
c
c******************************************************************************
	xini=0.005
	xfin=d
	convergencia=.false.
c
	i_ite=1
	do while((.not.convergencia).and.(i_ite.lt.100))
	  x=(xini+xfin)/2.
	  ec=x*ey/(d-x)
c       ContribuciÛn del concreto en la compresiÛn
	  if(op_diseno.eq."noconfinado") then
	      call f_comp_concret(ec,b,x,h,fc,eo,euc,
     s                             cc,mc,c1,c2)
c
	  else ! op_diseÒo = confinado
	      call f_comp_conconfi(ec,hc,x,h,fc,eo,rv,k,zm,fccmin,
     e                                  euc,eum,eccu,b,
     s                                  cc,mc)
	  end if
c
c        ContribuciÛn del acero m·s comprimido,asc, en la compresiÛn
	  eps=ec*(x-dc)/x
	  call esf_aceros(eps,fy,fsu,mesy,mesh,ey,esh,esm,
     s                       fps)
	  cs1=fps*asc
c
	  if (cs1.le.0.) then
	    f=0.
	  else if (cs1.gt.0.) then
		  if(op_diseno.eq."noconfinado")then
		    call esf_concret(eps,fc,eo,euc,
     s                          f)
		  else
			call esf_concret_confi(eps,fc,eo,k,zm,
     e                                  eccu,eum,
     s                                  f)
     		  end if
	  end if
c	  correcciÛn en la Fuerza de compresiÛn, a nivel del acero m·s comprimido
        aux1=asc*f
c
	  t=0.
	  aux2=0.
	  c_s=0.
c	Fuerzas a nivel de los aceros de las capas intermedias
	  do i=1,5
c	   Distancia al eje neutro de c/capa intermedia de acero(escalar)
	    d11(i)=x-di1(i)
c         Deformaciones unitarias al nivel de c/capa intermedia de acero
		es_interm(i)=d11(i)*ey/(d-x)
c		esfuerzos a nivel de los aceros intermedios
		call esf_aceros(es_interm(i),fy,fsu,mesy,mesh,ey,esh,esm,
     s                       fs_interm(i))
		if(dabs(fs_interm(i)).gt.fy) then
			fs_interm(i)=fy*es_interm(i)/dabs(es_interm(i))
		end if
c		Fuerza en Nivel de los aceros intermedios
		Fz_s_interm(i)=fs_interm(i)*as_interm(i)
c		C·lculo de la fuerza resultante tanto a compresiÛn como tracciÛn
c         a nivel de los aceros
		if(as_interm(i).eq.0.) then
			f_compre(i)=0.
		else if (as_interm(i).gt.0.) then
			if (es_interm(i).gt.0.) then
c				C·lculo del esfuerzo a compresiÛn en el nivel del acero
				if(op_diseno.eq."noconfinado")then
					call esf_concret(es_interm(i),fc,eo,euc,
     s                          f_compre(i))
				else if(op_diseno.eq."confinado")then
					call esf_concret_confi(es_interm(i),fc,eo,k,zm,
     e                                  eccu,eum,
     s                                  f_compre(i))
				end if
c				corecciÛn de la resultante de fuerza a compresiÛn
				aux2=aux2+as_interm(i)*f_compre(i)
c				Fuerza a compresiÛn a nivel de los aceros
				c_s=c_s+Fz_s_interm(i)
			else if (es_interm(i).le.0.) then
c				Fuerza a tracciÛn a nivel de los aceros
				t=t+Fz_s_interm(i)
				f_compre(i)=0.
			end if
c
		end if
	  end do
c
	  cs=c_s+cs1
	  comp=cc+cs-aux1-aux2
	  trac=-fy*ast+t
c
	  if(x.gt.d) then
		stop "eje neutro mayor que altura util"
c
	  else if(dabs(comp+trac).gt.(0.00001*dabs(trac))) then
			if(comp.lt.dabs(trac)) then
				xini=x
			else
				xfin=x
			end if
c
	  else
		convergencia=.true.
c
	  end if
	  i_ite=i_ite+1
c
	end do
	if (i_ite.eq.100)then
c           print*,'i =',i_ite
		stop "no hay convergencia en el calculo de mp bajo carga cero"
	end if
c
	if(dabs(eps).gt.ey)then
	 open(file="error",unit=05)
		write(05,*)"atencion: el acero en compresion esta"
		write(05,*)"fluyendo antes que el de traccion para"
		write(05,*)"carga axial cero"
		write(05,*)"el criterio utilizado para el calculo de Mp"
		write(05,*)"es la fluencia del acero en traccion"
	end if
c
	m_s=0.
	do i=1,5
	 m_s=m_s+(Fz_s_interm(i)-as_interm(i)*f_compre(i))*(h/2.-di1(i))
	end do
	mp=mc+cs1*(h/2.-dc)-aux1*(h/2.-dc)+ast*fy*(h/2.-dt)+m_s
	tabla_mp(2,1)=0.

      p=mp/1000.
	if (dabs(p).le.dabs(tabla_mcr(2,2))) then
	  tabla_mp(2,2)=1.02*tabla_mcr(2,2)
	else
        tabla_mp(2,2)=p
      endif

	tabla_mp(2,3)=ey*100./(d-x)/100.
c     print*,'mp2=',tabla_mp(2,2)
c	print*,'valor de x en mp2=',x
c	print*,'valor de curvatuta en mp2=',tabla_mp(2,3)
c	print*,'*********                     **********'
c
c******************************************************************************
c       TERCER VALOR DE LA TABLA: CALCULO DE LA CARGA Y MOMENTO PARA        
c	  PUNTO 3 MOMENTO PLASTICO
c       AMBOS ACEROS EN FLUENCIA Y SU CURVATURA CEDENTE
c******************************************************************************
c
	fs=-fy
c
	x=eo*d/(ey+eo)
c
	eps=(ey+eo)*(x-dc)/d
c
	ec=eo
c
	call esf_aceros(eps,fy,fsu,mesy,mesh,ey,esh,esm,fps)
c
c 
*****************************************************************************
c
c	C·lculo de la fuerza de compresiÛn en el concreto
	if(op_diseno.eq."noconfinado") then
		call f_comp_concret(ec,b,x,h,fc,eo,euc,
     s                             cc,mc,c1,c2)
	else
		call f_comp_conconfi(ec,hc,x,h,fc,eo,rv,k,zm,fccmin,
     e                                  euc,eum,eccu,b,
     s                                  cc,mc)
	end if
c     Esfuerzo en el concreto a nivel de la fibra de acero m·s comprimida
	if(op_diseno.eq."noconfinado")then
		call esf_concret(eps,fc,eo,euc,
     s                          f)
	else
		call esf_concret_confi(eps,fc,eo,k,zm,
     e                                  eccu,eum,
     s                                  f)
	end if
c     Fuerza de compresiÛn en el acero asc
	cs1=fps*asc
c
	aux2=0.
	m_s=0.
	c_s=0.
	t=0.
	do i=1,5
		d11(i)=x-di1(i)
		es_interm(i)=ey*d11(i)/(d-x)
		call esf_aceros(es_interm(i),fy,fsu,mesy,mesh,ey,esh,esm,
     s                       fs_interm(i))
c		Fuerza en Nivel de los aceros intermedios
		Fz_s_interm(i)=fs_interm(i)*as_interm(i)
c		C·lculo de la fuerza resultante tanto a compresiÛn como tracciÛn
c         a nivel de los aceros
		if (as_interm(i).eq.0.) then
			f_compre(i)=0.
		else if (as_interm(i).gt.0.) then
			if (es_interm(i).gt.0.) then
c				C·lculo del esfuerzo a compresiÛn en el nivel del acero
				if(op_diseno.eq."noconfinado")then
					call esf_concret(es_interm(i),fc,eo,euc,
     s                          f_compre(i))
				else if(op_diseno.eq."confinado")then
					call esf_concret_confi(es_interm(i),fc,eo,k,zm,
     e                                  eccu,eum,
     s                                  f_compre(i))
				end if
c				corecciÛn de la resultante de fuerza a compresiÛn

c				Fuerza a compresiÛn a nivel de los aceros
				c_s=c_s+Fz_s_interm(i)
			else if (es_interm(i).le.0.) then
c				Fuerza a tracciÛn a nivel de los aceros
				t=t+Fz_s_interm(i)
				f_compre(i)=0.
			end if
		end if
			aux2=aux2+as_interm(i)*f_compre(i)
			m_s=m_s+(Fz_s_interm(i)-as_interm(i)*
     1                       f_compre(i))*(h/2.-di1(i))
	end do
c
	cs=cs1+c_s
	comp=cc+cs-asc*f-aux2
	trac1=ast*fs
	trac=trac1+t
	pp=trac+comp
	mp=mc+cs1*(h/2.-dc)-trac1*(h/2.-dt)-asc*f*(h/2.-dc)+m_s
c
	tabla_mp(3,1)=-(pp/1000.)
c      print*,'Faxial para Mp3 =',tabla_mp(3,1)
	tabla_mp(3,2)=mp/1000.
c      print*,'Mp3 =',tabla_mp(3,2)
	tabla_mp(3,3)=ey*100./(d-x)/100.
c
c
c   cambio de los puntos de los diagramas de interacciÛn cuando la carga axial 
c   del punto Mp(3,1) es positiva
c
c
      if(tabla_mp(3,1).ge.0)then
	    Maux=tabla_mp(3,2) 
		Paux=-tabla_mp(3,1)
		tabla_mp(3,2)=tabla_mp(2,2)
		tabla_mp(3,1)=-tabla_mp(2,1)
		tabla_mp(2,2)=Maux
		tabla_mp(2,1)=Paux
	end if
c
c
c******************************************************************************
c       ULTIMO VALOR DE LA TABLA PARA CARGA AXIAL PURA Y                  
c	  PUNTO 4(ARRIBA ) MOMENTO PLASTICO
c       SU CURVATURA CEDENTE
c******************************************************************************
c
	if(op_diseno.eq."noconfinado")then
	   call esf_concret(ey,fc,eo,euc,
     s                          f)
	else
	   call esf_concret_confi(ey,fc,eo,k,zm,
     e                                  eccu,eum,
     s                                  f)
	end if
	Pplastico=((b*h-asc-ast-as_t_inter)*f+(ast+asc+as_t_inter)*fy)
c
      if (dabs(Pplastico).le.dabs(tabla_mcr(4,1)))then
	  tabla_mp(4,1)=-(1.02*tabla_mcr(4,1)/1000.)
	else
        tabla_mp(4,1)=-(Pplastico/1000.)
      endif
	tabla_mp(4,2)=0.
	tabla_mp(4,3)=0.
c
	return
	end

c
c******************************************************************************
c
c       CAL_MU
c       DIAGRAMA DE INTERACCION DEL DIAGRAMA DEL MOMENTO
c       ULTIMO
c       ESCRITA POR MARIA ELENA PERDOMO Y JULIO FLOREZ LOPEZ                
c       MODIFICADA POR RICARDO PIC”N, SCARLET Y JULIO FL”REZ 2021
c       CALCULO DE MOMENTO ULTIMO
c       INICIADA EL 21 DE JULIO DE 1997
c	  Modificada por: Nayive Jaramillo 2002
c******************************************************************************
c
	subroutine cal_mu(op_diseno,tabla_mp,
     e              h,fy,fsu,ey,esh,esm,
     e              mesy,mesh,ast,asc,w,as_interm,di1,d,dt,b,
     e              dc,rv,rh,ps,sh,fyh,k,zm,fc,eo,euc,eum,eccu,
     e              hc,fccmin,
     s              tabla_mu)
c
       IMPLICIT NONE
c******************************************************************************
c               DECLARACIONES LOCALES
c******************************************************************************
c
c      op_diseno=opciones de diseno para concreto confinad o no
c                 confinado
c       max_def=tamano de la matriz de deformaciones
c       esfuerzo_cero=esfuerzos al principio del paso
c       fr=modulo de rotura del concreto
c       inerc=inercia de la seccion
c       y=distancia del eje centroidal a la fibra mas traccionada
c       fy=esfuerzo de cedencia del acero de refuerzo
c       fsu=esfuerzo maximo del acero de refuerzo
c       ey=deformacion cedente del acero de refuerzo
c       esh=deformacion al final de la cedencia del acero
c       esm=deformacion maxima del acero de refuerzo
c       mesy=modulo de elasticidad del acero en la zona elastica
c       mesh=modulo 1eelasticidad del acero en la zona de endu-
c            recimiento.
c       ast=area de acero a traccion
c       asc=area de acero a compresion
c       d=altura util
c       dc=recubrimiento del acero a compresion
c       dt=recubrimiento del acero a traccion
c       y=distancia del eje centroidal a la fibra mas traccio-
c         nada
c       h=altura de la seccion
c       hc=ancho del nucleo confinado
c       rv=recubrimiento vertical de los estribos
c       k=factor de incremento de resistencia del concreto con-
c         finado.
c       zm=parametro que define la pendiente de la rama descen-
c          diente
c       fc=resistencia maxima a compresion del concreto
c       eo=deformacion del concreto para la resistencia maxima
c       euc=deformacion ultima del concreto no comfinado o defor-
c           macion al final de la rama descendente del concreto
c           confinado
c       eum=deformacion ultima del concreto confinado
c       l=longitud del elemento
c       n=fuerza axial
c       mcr=momento de agrietamiento
c       mp=momento plastico
c       mu=momento ultimo
c       fpu=rotula plastica ultima
c       x=profundidad del eje neutro
c       ec=deformacion del concreto en la fibra mas comprimida
c       eps=deformacion en la fibra por donde pasa la resultante
c           a compresion del acero de refuerzo
c       es=deformacion de la fibra mas traccionada
c       fs=esfuerzo del acero en la fibra mas traccionada
c       fps=esfuerzo del acero en la fibra por donde pasa la re-
c           sultante a compresion del acero de refuerzo
c       cs=resultante de la fuerza de compresion del acero de
c          refuerzo
c       cc=resultante de las fuerzas de compresion del concreto
c       mc=momento resistente debido a las fuerzas de compresion
c          del concreto
c       trac=resultante de las fuerzas de traccion
c       comp=resultante de las fuerzas de compresion
c       fiy=rotacion de la seccion correspondiente a la primera
c           fluencia del acero
c       fiu=rotacion de la seccion correspondiente a la rotura
c       lp=longitud de la articulacion plastica
c	  w=5= valor entero que indica en tamaÒo del vector aceros intermedios
c
c******************************************************************************
c
	integer w
	real*8 as_interm(w),di1(w),as_t_interm
	real*8 es_interm(5),d11(5),fs_interm(5),corec
	integer i
	real*8 m_s,fcomp_s(5),aux1(5),Fz_s_interm(5)
	real*8 t_s_interm,c_s_interm,cs1,trac1
	real*8 h,ast,asc,d,dc,dt
	real*8 hc,rv,rh,ps,sh,fyh,k,zm
	real*8 fy,fsu,ey,esh,esm,es,fs,eps,fps
	real*8 mesy,mesh
	real*8 fc,eo,euc,eum,ec,eccu
	real*8 x,x1,cc,mc,cs,comp,trac
	real*8 fccmin
	real*8 b,c1,c2,momen,momento
	character*14 op_diseno
	real*8 tabla_mp(4,3)
	real*8 tabla_mu(5,3)
	real*8 mu,curv_u,var
	real*8 m,curv
	real*8 f,pu,p,mom
	real*8 ec1
	real*8 maux
c
c******************************************************************************
c       PRIMER VALOR DE LA TABLA
c       FALLA DEL ACERO EN TRACCION                                       
c	  PTO. 1 (ABAJO ) MOMENTO ULTIMO
c******************************************************************************
c
	call esf_aceros(esm,fy,fsu,mesy,mesh,ey,esh,esm,
     s                         fps)
c
	as_t_interm=0.
	do i=1,5
		 as_t_interm=as_interm(i)+as_t_interm
	end do
      p=-(ast+asc+as_t_interm)*fps/1000.
	if (dabs(p).lt.dabs(tabla_mp(1,1))) then
	  tabla_mu(1,1)=(tabla_mp(1,1))
	else
        tabla_mu(1,1)=-p
      endif
	tabla_mu(1,2)=0.
	tabla_mu(1,3)=0.
c
c******************************************************************************
c       CALCULO DEL MOMENTO PARA CARGA AXIAL                             
c	  PTO.2 ( AXIAL = 0 ) MOMENTO ULTIMO
c       CERO Y SU CURVATURA
c******************************************************************************
	if (op_diseno.eq."noconfinado") then
		ec=euc
		var=0.
		call cal_mnc(ec,var,w,as_interm,di1,
     e                 euc,d,b,h,fc,eo,dc,dt,
     e                 fy,fsu,mesy,mesh,ey,esh,esm,asc,ast,
     s                 m,curv_u,x)
	else if (op_diseno.eq."confinado") then
			ec=eum
			var=0.
			call cal_mcc(ec,var,w,as_interm,di1,
     e                  euc,d,b,h,fc,eo,dc,dt,
     e                  eum,hc,rv,eccu,k,zm,
     e                  fy,fsu,mesy,mesh,ey,esh,esm,asc,ast,
     s                  m,curv_u,x,es)
c			call cal_m(ec,p,as_interm,di1,fccmin,
c     e                  euc,eccu,d,b,h,fc,eo,dc,dt,
c     e                  eum,hc,rv,k,zm,op_diseno,
c     e                  fy,fsu,mesy,mesh,ey,esh,esm,asc,ast,
c     s                  m,curv)
	end if
c
c
c
      if((m/1000.).le.(tabla_mp(2,2)))then
		tabla_mu(2,2)=tabla_mp(2,2)*1.02
c	print*,'cambi mu(p=0)=1.02mp(p=0)'
	else
		tabla_mu(2,2)=m/1000.
	end if
c      print*,'valor de x ultimo (carga= cero=)',x
      if(tabla_mp(2,1).gt.0)then
          tabla_mu(2,1)=tabla_mp(2,1)
      else
	    tabla_mu(2,1)=0.
	end if
	tabla_mu(2,3)=curv_u*100./100.
c******************************************************************************
c       CALCULO DEL MOMENTO PARA LA CARGA BALANCEADA DE FLUENCIA             
c           PTO.3
c******************************************************************************
c
      if (op_diseno.eq."noconfinado") then
		ec=euc
		var=-tabla_mp(3,1)*1000.
		call cal_mnc(ec,var,w,as_interm,di1,
     e                 euc,d,b,h,fc,eo,dc,dt,
     e                 fy,fsu,mesy,mesh,ey,esh,esm,asc,ast,
     s                 m,curv_u,x)
	else if (op_diseno.eq."confinado") then
			ec=eum
			var=-tabla_mp(3,1)*1000.
			call cal_mcc(ec,var,w,as_interm,di1,
     e                  euc,d,b,h,fc,eo,dc,dt,
     e                  eum,hc,rv,eccu,k,zm,
     e                  fy,fsu,mesy,mesh,ey,esh,esm,asc,ast,
     s                  m,curv_u,x,es)
	end if

c
	tabla_mu(3,1)=-var/1000.
      momen=m/1000.
	if (momen.le.tabla_mp(3,2)) then
          tabla_mu(3,2)=1.02*tabla_mp(3,2)
c		print*,'cambio en el mu para la carga balanceada de MU'
c		print*,'tabla_mu(3,2)',tabla_mu(3,2)
c		print*,'Mu calculado',momen
c		print*,'tabla_mp(3,2)',tabla_mp(3,2)
	else
	    tabla_mu(3,2)=momen
	endif
c
c******************************************************************************
c       CALCULO DEL MOMENTO ULTIMO BALANCEADO, SU CARGA AXIAL Y SU CURVATURA 
c           PTO. 4
c******************************************************************************
c
	if(op_diseno.eq."noconfinado")then
		x=euc*d/(euc+ey)
		ec=euc
		eps=euc*(x-dc)/x
		call f_comp_concret(ec,b,x,h,fc,eo,euc,
     s                             cc,mc,c1,c2)
		curv_u=euc/x
	else if(op_diseno.eq."confinado")then
		x1=eum*(d-rv)/(eum+ey)
c		x1=eum*(d-rv)/(eum+0.11)
c		ec=eum*(x1+rv)/x1
		ec=eum
		eps=eum*(x1+rv-dc)/x1
		x=x1+rv
		call f_comp_conconfi_mu(ec,hc,x1,h,fc,eo,rv,k,
     e                               euc,eum,eccu,b,
     s                               cc,mc)
		curv_u=eum/x1
c
	end if
	call esf_aceros(eps,fy,fsu,mesy,mesh,ey,esh,esm,
     s                         fps)

	trac1=-ast*fy
	cs1=asc*fps
	if(op_diseno.eq."noconfinado")then
	   call esf_concret(eps,fc,eo,euc,
     s                          f)
	else
	   call esf_concret_confi(eps,fc,eo,k,zm,
     e                                  eccu,eum,
     s                                  f)
	end if
c
	m_s=0.
	t_s_interm=0.
	c_s_interm=0.
	corec=0.
	do i=1,5
		d11(i)=x-di1(i)
		if(op_diseno.eq."noconfinado") then
			es_interm(i)=(ey+ec)*d11(i)/d
			if (es_interm(i).gt.0.) then
				call esf_concret(es_interm(i),fc,eo,euc,
     s                          fcomp_s(i))
			else if (es_interm(i).le.0.) then
				fcomp_s(i)=0.
			end if
		else if (op_diseno.eq."confinado") then
			es_interm(i)=(ey+ec)*d11(i)/(d-rv)
			if (es_interm(i).gt.0.) then
				call esf_concret_confi(es_interm(i),fc,eo,k,zm,
     e                                  eccu,eum,
     s                                  fcomp_s(i))
			else if (es_interm(i).le.0.) then
				fcomp_s(i)=0.
			end if
		end if
c
		aux1(i)=as_interm(i)*fcomp_s(i)
		corec=corec+aux1(i)
c
		call esf_aceros(es_interm(i),fy,fsu,mesy,mesh,ey,esh,esm,
     s                         fs_interm(i))
c
c		if(dabs(es_interm(i)).gt.esm) then
c			fs_interm(i)=fsu*dabs(es_interm(i))/(es_interm(i))
c		end if
c
		Fz_s_interm(i)=as_interm(i)*fs_interm(i)
		m_s=m_s+(Fz_s_interm(i)-as_interm(i)*fcomp_s(i))*(h/2.-di1(i))
		if (Fz_s_interm(i).gt.0.) then
			c_s_interm=c_s_interm+Fz_s_interm(i)
		else if (Fz_s_interm(i).le.0.) then
			t_s_interm=t_s_interm+Fz_s_interm(i)
		end if
	end do
	cs=cs1+c_s_interm
	comp=cc+cs-asc*f-corec
	trac=trac1+t_s_interm
	pu=trac+comp
	mu=mc+cs1*(h/2.-dc)-trac1*(h/2.-dt)-asc*f*(h/2.-dc)+m_s
c
	tabla_mu(4,1)=-pu/1000.
c	print*,'carga axial 4,1=',tabla_mu(4,1)
      momento=mu/1000.
c	print*,'CALCULADO=',momento
c	print*,'tabla_mp(3,2)=',tabla_mp(3,2)
      if (momento.le.tabla_mp(3,2)) then
          tabla_mu(4,2)=1.03*tabla_mp(3,2)
	else
	    tabla_mu(4,2)=momento
	end if
      tabla_mu(3,3)=tabla_mu(2,3)/7.
	tabla_mu(4,3)=curv_u*100./100.
c
c******************************************************************************
c       ULTIMO VALOR DE LA TABLA                                           
c	  PTO. 5 ( ARRIBA ) MOMENTO ULTIMO
c******************************************************************************
c
c       CALCULO POR FALLA DEL CONCRETO EN COMPRESION
c
	if(op_diseno.eq."noconfinado")then
		f=0.85*fc
		es=eo
		call esf_aceros(es,fy,fsu,mesy,mesh,ey,esh,esm,
     s                         fs)
	    pu=((b*h-asc-ast-as_t_interm)*f+
     1		(ast+asc+as_t_interm)*fs)/1000.
	else
		f=k*fc
		es=k*eo
		call esf_aceros(es,fy,fsu,mesy,mesh,ey,esh,esm,
     s                         fs)
		pu=((b*h-asc-ast-as_t_interm)*f+
     1			(ast+asc+as_t_interm)*fs)/1000.
c
	end if
	if (pu.lt.-tabla_mp(4,1)) then
	   tabla_mu(5,1)=1.02*tabla_mp(4,1)
	else
	   tabla_mu(5,1)=-pu
      endif
	tabla_mu(5,2)=0.
	tabla_mu(5,3)=0.
c
c       print*,'Pplastico=',tabla_mp(4,1)
	return
	end
c******************************************************************************
c
c       SUBRUTINA ESF_CONCRET PARA CUALQUIER MODELO
c       CALCULO DEL ESFUERZO DE COMPRESION EN EL CONCRETO
c       NO CONFINADO CONSIDERANDO EL DIAGRAMA ESFUERZO-DEFORMACION
c       PROPUESTO POR HOGNESTAD.
c
c
c
c******************************************************************************
c
	subroutine esf_concret(ec,fc,eo,euc,
     s                          f)
c******************************************************************************
c               DECLARACIONES LOCALES
c******************************************************************************
c
       IMPLICIT NONE
c       ec=deformacion a compresion en cualquier fibra
c       fc=resistencia maxima a compresion del concreto
c       eo=deformacion del concreto para la resistencia maxima
c       euc=deformacion maxima del concreto no confinado
c
c******************************************************************************
	real*8 ec,eo,fc,euc
	real*8 f
c
c******************************************************************************
c  Calculo del esfuerzo a compresion para una deformacion cualquiera
c******************************************************************************
c				DIAGRAMA ESFUERZO-DEFORMACION
c				PROPUESTO POR HOGNESTAD.

	 if(ec.lt.eo)then
		f=0.85*fc*(2.*ec/eo-(ec/eo)**2)
	  else if(ec.lt.euc)then
		f=0.85*fc*(1.-0.15*(ec-eo)/(euc-eo))
	  else
		f=0.
	  end if
c
	return
	end
c******************************************************************************
c
c       SUBRUTINA ESF_CONCRET_CONFI PARA CUALQUIER MODELO
c       CALCULO DEL ESFUERZO DE COMPRESION EN EL CONCRETO
c       CONFINADO CONSIDERANDO EL DIAGRAMA ESFUERZO-DEFORMACION
c       PROPUESTO POR KENT Y PARK.
c
c******************************************************************************
c
	subroutine esf_concret_confi(ec,fc,eo,k,zm,
     e                                  eccu,eum,
     s                                  f)
c******************************************************************************
c               DECLARACIONES LOCALES
c******************************************************************************
c
       IMPLICIT NONE
c       ec=deformacion a compresion en cualquier fibra
c       fc=resistencia maxima a compresion del concreto
c       eo=deformacion del concreto para la resistencia maxima
c       euc=deformacion maxima del concreto no confinado
c
c******************************************************************************
	real*8 ec,eo,fc,eccu,eum
	real*8 k,zm
	real*8 f
c
c******************************************************************************
c  Calculo del esfuerzo a compresion para una deformacion cualquiera
c******************************************************************************
c			DIAGRAMA ESFUERZO-DEFORMACION
c			PROPUESTO POR KENT Y PARK.
	if(ec.le.k*0.002)then
		f=k*fc*(2.*ec/(k*0.002)-(ec/(k*0.002))**2)
	  else if(ec.le.eccu)then
		f=k*fc*(1.-zm*(ec-k*0.002))
	  else if(ec.lt.eum)then
		f=0.2*k*fc
	  else
		f=0.
c
	end if
c
	return
	end
c
c******************************************************************************
c
c       SUBRUTINA F_COMP_CONCRET PARA CUALQUIER MODELO
c       CALCULO DE LA FUERZA DE COMPRESION RESULTANTE EN EL CONCRE-
c       TO PARA SECCIONES RECTANGULARES DE CONCRETO ARMADO NO CON-
c       FINADO CONSIDERANDO EL DIAGRAMA ESFUERZO-DEFORMACION PRO-
c       PUESTO POR HOGNESTAD.
c
c******************************************************************************
c
	subroutine f_comp_concret(ec,b,x,h,fc,eo,euc,
     s                            cc,mc,c1,c2)
c
       IMPLICIT NONE
c******************************************************************************
c               DECLARACIONES LOCALES
c******************************************************************************
c
c       ec=deformacion del concreto en la fibra mas comprimida
c       b=ancho de la seccion
c       x=profundidad del eje neutro
c       d=altura util
c       fc=resistencia maxima a compresion del concreto
c       eo=deformacion del concreto para la resistencia maxima
c       euc=deformacion maxima del concreto no confinado
c       cc=resultante de las fuerzas de compresion del concreto
c       mc=momento resistente debido a las fuerzas de compresion
c          del concreto
c
c******************************************************************************
c
	real*8 ec,eo,fc,euc,ya
	real*8 b,x,h
	real*8 cc,mc
	real*8 c1,c2,y1,y2,f,auxiliar
c
c******************************************************************************
c       Calculo de la fuerza resultante de compresion del concreto    
*********
c******************************************************************************
c
	if(ec.le.eo) then
	   c1=0.85*fc*b*x*(3.*ec*eo-(ec*ec))/(3.*(eo*eo))
	   y1=x*(8.*eo-3.*ec)/(4.*(3.*eo-ec))
	   c2=0
	   y2=0
c
	else if((ec.gt.eo).and.(ec.le.euc)) then
	   c1=(2.*0.85*fc*b*x*eo)/(3.*ec)
	   y1=(x*5.*eo)/(8.*ec)
	   auxiliar=0.85*fc*(1.-0.15*(ec-eo)/(euc-eo))
	   c2=x*b*(auxiliar+0.85*fc)*(ec-eo)/(2.*ec)
	   ya=(3.*(auxiliar+0.85*fc))
	   y2=x*((eo/ec)+(1.-(eo/ec))*(2.*0.85*fc+auxiliar)/ya)
c
	end if
c
	      cc=c1+c2
c******************************************************************************
c       Calculo del momento resistente debido a las fuerzas de com-
c       presion del concreto.
c******************************************************************************
c
	      mc=c1*(h/2.-x+y1)+c2*(h/2.-x+y2)
c
	return
	end
c******************************************************************************
c
c       SUBRUTINA F_COMP_CONCONFI
c
c       CALCULO DE LA FUERZA RESULTANTE DE COMPRESION EN EL CONCRE-
c       TO PARA SECCIONES RECTANGULARES DE CONCRETO CONFINADO CON-
c       SIDERANDO EL DIAGRAMA ESFUERZO-DEFORMACION PROPUESTO POR
c       KENT Y PARK MODIFICADO, se considera que la secciÛn tiene
c	  una zona confinada y otra no confinada (se resuelve por
c       HOGNESTAD)
C
C
c******************************************************************************
C
	subroutine f_comp_conconfi(ec,hc,x,h,fc,eo,rv,k,zm,fccmin,
     e                                 euc,eum,eccu,b,
     s                                 cc,mc)
c
       IMPLICIT NONE
c******************************************************************************
c               DECLARACIONES LOCALES
c******************************************************************************
c
c       hc=ancho del nucleo confinado
c       x=profundidad del eje neutro
c       fc=resistencia maxima a compresion del concreto
c       eo=deformacion del concreto para la resistencia maxima
c       rv=recubrimiento vertical de los estribos
c       k=factor de incremento de resistencia del concreto confi-
c         nado
c       zm=parametro que define la pendiente de la rama descen-
c          dente del concreto confinado
c       euc=deformacion al final de la rama descendente del concre-
c           to  no confinado
c       eum=deformacion maxima del concreto confinado
c       cc=resultante de las fuerzas de compresion del concreto
c       mc=momento resistente debido a las fuerzas de compresion
c          del concreto
c       eccu=deformaciÛn ultima del concreto confinado
c
c******************************************************************************
c
	real*8 ec,eo,euc,eum,e1,eccu,ecc,ecco,eh1,emax1
	real*8 hc,x,h,yo,x1,b,a,b1
	real*8 fccmin,fcc,fc,fy,f,fh,k,zm,yu,w
	real*8 rv,s,t,sh,Yn,a1,a2,q
	real*8 mc,m1a,m1b,Mcnc
	real*8 comp1,comp2,comp3,brazo1,brazo2,brazo3
      real*8 cc,c1a,c1b,Ccnc,c1aa,c2aa,c1bb,c2bb
c
c******************************************************************************
c    Calculo de la fuerza resultante de compresion del concreto de la parte 
c    NO CONFINADA
c******************************************************************************
c     HAY QUE REVISAR SI SOBRAN ALGUNAS VARIABLES DECLARADAS
      ecc=ec*(x-rv)/x
	x1=x-rv
c
	if (ec.le.eo) then
	    call f_comp_concret(ec,b,x,h,fc,eo,euc,
     s                             c1a,m1a,c1aa,c2aa)
          call f_comp_concret(ecc,hc,x1,h,fc,eo,euc,
     s                             c1b,m1b,c1bb,c2bb)
	    Ccnc=c1a-c1b
          Mcnc=m1a-m1b
c
      else if ((ec.gt.eo).and.(ec.le.euc)) then
	        call f_comp_concret(ec,b,x,h,fc,eo,euc,
     s                             c1a,m1a,c1aa,c2aa)
              call f_comp_concret(ecc,hc,x1,h,fc,eo,euc,
     s                             c1b,m1b,c1bb,c2bb)
	        Ccnc=c1a-c1b
              Mcnc=m1a-m1b
c
           else if ((ec.gt.euc).and.(ecc.lt.euc)) then
	           Yn=x*euc/ec
                 call f_comp_concret(euc,b,Yn,h,fc,eo,euc,
     s                             c1a,m1a,c1aa,c2aa)
                 call f_comp_concret(ecc,hc,x1,h,fc,eo,euc,
     s                             c1b,m1b,c1bb,c2bb)

                 Ccnc=c1a-c1b
                 Mcnc=m1a-m1b
	           else if ((ec.gt.euc).and.(ecc.gt.euc)) then
	               Yn=x*euc/ec
                     call f_comp_concret(euc,b,Yn,h,fc,eo,euc,
     s                             c1a,m1a,c1aa,c2aa)
                     call f_comp_concret(euc,hc,Yn,h,fc,eo,euc,
     s                             c1b,m1b,c1bb,c2bb)
                     q=(ec-euc)/euc
                     Ccnc=c1a-c1b
                     Mcnc=(m1a-(c1aa+c2aa)*q*Yn)-(m1b-(c1bb+c2bb)*q*Yn)
	end if
c******************************************************************************
c       Calculo de la fuerza resultante de compresion del concreto de la 
c	  parte  CONFINADA (Kent y Park)
c******************************************************************************
	ecco=0.002*k
	if (ecc.le.ecco) then
	    comp1=(k*fc)*hc*(x-rv)*(3.*ecco*ecc-ecc*ecc)/(3.*ecco*ecco)
	    brazo1=(x-rv)*(8.-3.*ecc/ecco)/(4.*(3.-ecc/ecco))
          comp2=0
          brazo2=0
	    comp3=0
          brazo3=0
	else if ((ecc.gt.ecco).and.(ecc.lt.eccu)) then
c
	        fcc=(k*fc)-(((k*fc)-fccmin)/(eccu-ecco))*(ecc-ecco)
              comp1=2.*k*fc*hc*(x-rv)*ecco/(3.*ecc)
	        brazo1=5.*ecco*(x-rv)/(8.*ecc)
	        comp2=(hc*(x-rv)/(2*ecc))*(fcc+k*fc)*(ecc-ecco)
	        brazo2=(x-rv)*((ecco/ecc)+(((1.-ecco/ecc)*(2.*k*fc)+
     1                  fcc))/(3.*(fcc+k*fc)))
	        comp3=0
	        brazo3=0
	     else if ((ecc.ge.eccu).and.(ecc.le.eum)) then
                  comp1=2.*k*fc*hc*(x-rv)*(ecco/ecc)/3.
	            brazo1=5.*ecco*(x-rv)/(8.*ecc)
	            comp2=hc*(x-rv)*(1.2*k*fc)*(eccu-ecco)/(2*eccu)
	            brazo2=(x-rv)*((ecco/eccu)+((1.-ecco/eccu)*
     1                     (2.2*k*fc))/(3.*1.2*k*fc))
	            yu=(x-rv)*eccu/ecc
	            comp3=0.2*k*fc*hc*(x-rv-yu)
              	brazo3=0.5*(x-rv+yu)
	    endif
c******************************************************************************
c       Calculo del momento resistente debido a las fuerzas de compresion 
c
c	  del  concreto
c******************************************************************************
c
	mc=comp1*(h/2.-x+brazo1)+comp2*(h/2.-x+brazo2)+comp3*
     1             (h/2.-x+brazo3)+Mcnc
	cc=comp1+comp2+comp3+Ccnc
c
	return
	end
c******************************************************************************
c
c       SUBRUTINA ESF_ACEROS PARA CUALQUIER MODELO
c       CALCULO DE LOS ESFUERZOS EN LOS ACEROS A TRACCION O COMPRE-
c       SION A PARTIR DE LAS DEFORMACIONES EN LOS CENTROS DE GRAVE-
c       DAD DE LAS AREAS RESPECTIVAS DE ACERO, CONSIDERANDO DIAGRA-
c       MA ESFUERZO-DEFORMACION CON ZONA DE ENDURECIMIENTO.
c	Fsu !           * * **
c		!	      *       *
c		!       *         *
c      Fy !    ***		  *
c		!   *			  *
c		!  *              *
c		! *               *
c		!*                *
c         !-----------------***********
c			ey	esh		esm
c
c
c******************************************************************************
c
	subroutine esf_aceros(es,fy,fsu,mesy,mesh,ey,esh,esm,
     s                     fs)
c
       IMPLICIT NONE
c******************************************************************************
c               DECLARACIONES LOCALES
c******************************************************************************
c
c       es=deformacion en el centro de gravedad del area de acero
c          a traccion o compresion
c       fy=esfuerzo de cedencia del acero de refuerzo
c       fsu=esfuerzo maximo del acero de refuerzo
c       mesy=modulo de elasticidad del acero en la zona elastica
c       mesh=modulo de elasticidad del acero en la zona de endure-
c            cimiento
c       ey=deformacion en la cedencia del acero de refuerzo
c       esh=deformacion al final de la cedencia del acero de re-
c           fuerzo
c       esm=deformacion del acero para el esfuerzo maximo
c       fs=esfuerzo en el centro de gravedad del area de acero a
c          traccion o a compresion
c
c******************************************************************************
c
	real*8 es,ey,esh,esm,fy,fsu,fs
	real*8 mesy,mesh,niu
c
c******************************************************************************
c
	if(dabs(es).le.ey) then
	      fs=es*mesy
	else if((dabs(es).gt.ey).and.(dabs(es).le.esh)) then
	      fs=fy*es/dabs(es)
	else if((dabs(es).le.esm)) then
	      niu=(dabs(es)-esh)/(esm-esh)
	      fs=(fy+(fsu-fy)*(2.*niu-(niu*niu)))*es/dabs(es)
	else
		  fs=0.
	end if
c
	return
	end
c
c******************************************************************
c
c
c       CALCULO DEL DIAGRAMA DE INTERACCION DE CURVATURAS PLASTICAS
c       ULTIMAS (fpu=fu-fy)                                                
c       CALCULO DE ROTACI”N ULTIMA
c
c**********************************************************************
c
	subroutine cal_x(tabla_mp,tabla_mu,
     s                  tabla_x)
c
       IMPLICIT NONE
c********************************************************************
c
	real*8 tabla_mp(4,4),tabla_mu(5,3),tabla_x(4,2),var,op
c
c********************************************************************
c
      op=10.
	tabla_x(1,1)=tabla_mu(1,1)
	tabla_x(1,2)=0.
c
	tabla_x(2,1)=0
	if ((tabla_mu(2,3)/op).gt.(tabla_mp(2,3))) then
	   tabla_x(2,2)=tabla_mu(2,3)/op-tabla_mp(2,3)
	else
         tabla_x(2,2)=tabla_mu(2,3)-tabla_mp(2,3)
	end if
c
c
	tabla_x(3,1)=tabla_mu(4,1)
	call interpolacion_x(tabla_mu(4,1),4,tabla_mp,4,
     s                          var)
c
      if ((tabla_mu(4,3)/op).gt.(var)) then
	   tabla_x(3,2)=tabla_mu(4,3)/op-var
	else
	   tabla_x(3,2)=tabla_mu(4,3)-var
	end if
c
	tabla_x(4,1)=tabla_mu(5,1)
	tabla_x(4,2)=0.
c
	return
	end
c
c**********************************************************************
c       CALCULO DEL MOMENTO Y LA CURVATURA PARA UNA FUERZA AXIAL Y UNA
c       DEFORMACION DEL CONCRETO DADA CUANDO EL CONCRETO SE CONSIDERA
c       NO CONFINADO (ec=deformacion en la fibra m·s comprimida,euc)
c       Elaborada por: Nayive Jaramillo
c
c**********************************************************************
c
	subroutine cal_mnc(ec,p,w,as_interm,di1,
     e                   euc,d,b,h,fc,eo,dc,dt,
     e                   fy,fsu,mesy,mesh,ey,esh,esm,asc,ast,
     s                   m,curv,x)
c
       IMPLICIT NONE
c**********************************************************************
c       DECLARACIONES
c**********************************************************************
c
c       ec=deformacion en el concreto
c       p=fuerza axial sobre la seccion
c       euc,d,b,h,fc,eo,dc,dt=datos para concreto no confinado
c       fy,fsu,mesy,mesh,ey,esh,esm,asc,ast=datos para el acero
c       m=momento sobre la seccion
c       curv=curvatura
c       fs_interm(i)=esfuezo a nivel de los aceros intermedios
c       f_comp_si(5)=esfuerzo en el concreto a nivel del acero intermedios
c	  (si el ‡cero est· sometido a compresiÛn ----> f_comp_si(i) >0.;
c	   si el acero se encuantra traccionado   ----> f_comp_si(i)=0.)
c
c**********************************************************************
c
	integer w
	real*8 as_interm(w),di1(w)
	real*8 ec,p
	real*8 euc,d,b,h,fc,eo,dc,dt
	real*8 fy,fsu,mesy,mesh,ey,esh,esm,asc,ast
	real*8 m,curv
c
c**********************************************************************
c       LOCALES
c**********************************************************************
c
	real*8 xini,xfin,x,es,eps,cs,cc,f,comp,trac,fs,fps,mc,c1,c2
	logical convergencia
	integer i_ite,i
	real*8 f_comp_si(5),Fz_s_interm(5),d11(5)
	real*8 aux1,correc,m_s,t_s_interm,c_s_interm
	real*8 cs1,es_interm(5),fs_interm(5),t,trac1
c
c**********************************************************************
c
	xini=0.005
	xfin=d
	convergencia=.false.
c
	i_ite=1
	do while((.not.convergencia).and.(i_ite.lt.50))
	  x=(xini+xfin)/2.
c
	  call f_comp_concret(ec,b,x,h,fc,eo,euc,
     s                             cc,mc,c1,c2)
c
	  eps=ec*(x-dc)/x
	  call esf_aceros(eps,fy,fsu,mesy,mesh,ey,esh,esm,
     s                       fps)
	  cs1=fps*asc
c
	  if (cs1.le.0.) then
		f=0.
	  else if (cs1.gt.0.) then
		  call esf_concret(eps,fc,eo,euc,
     s                          f)
	  end if
c
	  aux1=asc*f
	  correc=0.
	  m_s=0.
	  c_s_interm=0.
	  t=0.
	  do i=1,5
		d11(i)=x-di1(i)
		es_interm(i)=ec*d11(i)/x
		call esf_aceros(es_interm(i),fy,fsu,mesy,mesh,ey,esh,esm,
     s                       fs_interm(i))
c		Fuerza en Nivel de los aceros intermedios
		Fz_s_interm(i)=fs_interm(i)*as_interm(i)
c		C·lculo de la fuerza resultante tanto a compresiÛn como tracciÛn
c         a nivel de los aceros
		if (es_interm(i).gt.0.) then
c			C·lculo del esfuerzo a compresiÛn en el nivel del acero
			call esf_concret(es_interm(i),fc,eo,euc,
     s                          f_comp_si(i))
c			corecciÛn de la resultante de fuerza a compresiÛn
			correc=correc+as_interm(i)*f_comp_si(i)

c			Fuerza a compresiÛn a nivel de los aceros
			c_s_interm=c_s_interm+Fz_s_interm(i)
		else if (es_interm(i).le.0.) then
c			Fuerza a tracciÛn a nivel de los aceros
	        t=t+Fz_s_interm(i)
			f_comp_si(i)=0.
		end if
	end do
c
	  cs=cs1+c_s_interm
	  comp=cc+cs-aux1-correc
c
	  es=ec*(d-x)/x
	  call esf_aceros(es,fy,fsu,mesy,mesh,ey,esh,esm,
     s                       fs)
	  if (es.gt.esm) then
		fs=fsu
	  end if
	  trac1=-fs*ast
c
	  trac=trac1+t
	  if(x.gt.d) then
		stop "eje neutro mayor que altura util caso  no confinado"

	  else if(dabs(comp+trac-p).gt.(0.000000001*dabs(trac))) then
		if(comp.lt.dabs(trac-p)) then
c		   Aumenta la profundidad del Eje Neutro.
		   xini=x
		else
c		   Disminuye la profundidad del Eje Neutro.
		   xfin=x
		end if
	  else
		convergencia=.true.
	  end if
	  i_ite=i_ite+1
	end do
	if(i_ite.eq.100)then
		stop "no hay convergencia bajo carga cero en mu,(no confinado)"
	end if
c
	do i=i,5
		m_s=m_s+(Fz_s_interm(i)-as_interm(i)*f_comp_si(i))*
     1               (h/2.-di1(i))
	end do
	m=mc+cs1*(h/2.-dc)-trac1*(h/2.-dt)-asc*f*(h/2.-dc)+m_s
	curv=ec/x
c
	return
	end
c**********************************************************************
c
c      CALCULO DEL MOMENTO Y LA CURVATURA PARA UNA FUERZA AXIAL Y UNA
c      DEFORMACION DEL CONCRETO DADA en el caso que el concreto se
c      CONSIDERA CONFINDADO (ec=deformacion en la fibra mas conprimida,eum,
c	 a nivel del nucleo confinado )
c      Nayive Jaramillo S.
c**********************************************************************
c
	subroutine cal_mcc(ec,p,w,as_interm,di1,
     e                  euc,d,b,h,fc,eo,dc,dt,
     e                  eum,hc,rv,eccu,k,zm,
     e                  fy,fsu,mesy,mesh,ey,esh,esm,asc,ast,
     s                  m,curv,x,e1)
C
       IMPLICIT NONE
c**********************************************************************
c       DECLARACIONES
c**********************************************************************
C
c       ec=deformacion en el concreto
c       p=fuerza axial sobre la seccion
c       euc,d,b,h,fc,eo,dc,dt=datos para concreto no confinado
c       eum,hc,rv,k,zm,=datos para concreto confinado
c       fy,fsu,mesy,mesh,ey,esh,esm,asc,ast=datos para el acero
c       m=momento sobre la seccion
c       curv=curvatura
c
c       x=distancia desde el eje neutro hasta la fibra de concreto m·s
c	  comprimida del nucleo confinado
c
c	  e1= deformaciÛn el acero m·s traccionado
c**********************************************************************
c
	integer w
	real*8 as_interm(w),di1(w)
	real*8 ec,p
	real*8 euc,d,b,h,fc,eo,dc,dt
	real*8 eum,hc,rv,eccu,k,zm
	real*8 fy,fsu,mesy,mesh,ey,esh,esm,asc,ast
	real*8 m,curv,x,e1
c
c**********************************************************************
c       LOCALES
c**********************************************************************
c
	real*8 x1
	real*8 xini,xfin,es,eps,cs,cc,f,comp,trac,fs,fps,mc
	logical convergencia
	integer i_ite,i
	real*8 c1a,m1a,c1aa,c2aa,c1b,m1b,c1bb,c2bb,Yn,q,Ccnc,Mcnc
	real*8 ecc,ecco
	real*8 yu
	real*8 comp1,comp2,comp3
	real*8 brazo1,brazo2,brazo3
	real*8 cs1,t,m_s,trac1
	real*8 aux1,aux2,c_s,f_s(5),f_c_s(5)
	real*8 es_interm(5),Fz_s_interm(5),d11(5)
c
c**********************************************************************
c
	xini=0.005
	xfin=d
	convergencia=.false.
c
	i_ite=1
	do while((.not.convergencia).and.(i_ite.lt.100))
		x1=(xini+xfin)/2.
c		ContribuciÛn del concreto en la fuerza a compresiÛn
		call f_comp_conconfi_mu(ec,hc,x1,h,fc,eo,rv,k,
     e                                 euc,eum,eccu,b,
     s                                 cc,mc)
c
		eps=ec*(x1+rv-dc)/x1
		call esf_aceros(eps,fy,fsu,mesy,mesh,ey,esh,esm,
     s                       fps)
     	    cs1=fps*asc
c
		if (cs1.gt.0.) then
			call esf_concret_confi(eps,fc,eo,k,zm,
     e                                  eccu,eum,
     s                                  f)
		else if (cs1.le.0.) then
			f=0.
		end if
		aux1=asc*f
		aux2=0.
		c_s=0.
		t=0.
	    m_s=0.
c		Fuerzas a nivel de los aceros de las capas intermedias
		do i=1,5
c			Distancia al eje neutro de c/capa intermedia de
c			acero(escalar + o -)
c
			d11(i)=x1+rv-di1(i)
c
c             Deformaciones unitarias al nivel de c/capa
c			intermedia de acero
			es_interm(i)=d11(i)*ec/x1
c		    esfuerzos a nivel de los aceros intermedios
			call esf_aceros(es_interm(i),fy,fsu,mesy,mesh,ey,esh,esm,
     s                       f_s(i))
c**********************************************************
c			if((dabs(es_interm(i)).gt.esm).and.
c      1                                  (dabs(f_s(i)).le.fsu)) then
c				f_s(i)=fsu*es_interm(i)/dabs(es_interm(i))
c			end if
c*************************************************************
c		    Fuerza en Nivel de los aceros intermedios
			Fz_s_interm(i)=f_s(i)*as_interm(i)
c		    C·lculo de la fuerza resultante tanto a compresiÛn como tracciÛn
c             a nivel de los aceros
			if (es_interm(i).gt.0.) then
c			    C·lculo del esfuerzo a compresiÛn en el nivel del acero
				call esf_concret_confi(es_interm(i),fc,eo,k,zm,
     e                                  eccu,eum,
     s                                  f_c_s(i))
c			    Fuerza a compresiÛn a nivel de los aceros
				c_s=c_s+Fz_s_interm(i)
c        	        corecciÛn de la resultante de fuerza a compresiÛn
				aux2=aux2+as_interm(i)*f_c_s(i)
			else if (es_interm(i).le.0.) then
c			    Fuerza a tracciÛn a nivel de los aceros
				t=t+Fz_s_interm(i)
				f_c_s(i)=0.
			end if
		end do
c
		cs=cs1+c_s
		comp=cc+cs-aux1-aux2
c
		es=ec*(d-rv-x1)/x1
		call esf_aceros(es,fy,fsu,mesy,mesh,ey,esh,esm,
     s                       fs)
		if (es.gt.esm) then
			es=esm
			call esf_aceros(es,fy,fsu,mesy,mesh,ey,esh,esm,
     s                       fs)
		end if
		trac1=-fs*ast
c
		trac=trac1+t
		x=x1+rv
		if((x).gt.d) then
			stop "eje neutro mayor que altura util, caso confinado"
c
		else if(dabs(comp+trac-p).gt.(0.0000001*dabs(trac))) then
			if(comp.lt.dabs(trac-p)) then
				xini=x1
			else
				xfin=x1
			end if

	  else
c
		convergencia=.true.
c
	  end if
	  i_ite=i_ite+1
	end do
	if(i_ite.eq.100)then
		stop "no hay convergencia en el calculo de mu bajo carga cero"
	end if
c
	do i=1,5
		m_s=m_s+(Fz_s_interm(i)-as_interm(i)*f_c_s(i))*(h/2-di1(i))
	end do
	m=mc+cs1*(h/2.-dc)-trac1*(h/2.-dt)-asc*f*(h/2.-dc)+m_s
	curv=ec/x1
c
      e1=es
	return
	end
c******************************************************************************
c
c       SUBRUTINA F_COMP_CONCONFI_mu
c
c       CALCULO DE LA FUERZA RESULTANTE DE COMPRESION EN EL CONCRE-
c       TO PARA SECCIONES RECTANGULARES DE CONCRETO CONFINADO, en el
c       diagrama de momento ultimo, CONSIDERANDO EL DIAGRAMA
c       ESFUERZO-DEFORMACION PROPUESTO POR KENT Y PARK MODIFICADO,
c       se considera que la secciÛn tiene una zona confinada y otra
c	  no confinada (se resuelve por HOGNESTAD), ecc=ec=eum
c
c      Nayive Jaramillo S.
c******************************************************************************
C
	subroutine f_comp_conconfi_mu(ec,hc,x,h,fc,eo,rv,k,
     e                                 euc,eum,eccu,b,
     s                                 cc1,mc1)
c
       IMPLICIT NONE
c******************************************************************************
c               DECLARACIONES LOCALES
c******************************************************************************
c       hc=ancho del nucleo confinado
c
c       x=distancia desde el eje neutro hasta la fibra m·s comprimida del 
c		nucleo confinado
c
c       fc=resistencia maxima a compresion del concreto
c       eo=deformacion del concreto para la resistencia maxima
c       rv=recubrimiento vertical de los estribos
c       k=factor de incremento de resistencia del concreto confi-
c         nado
c       zm=parametro que define la pendiente de la rama descen-
c          dente del concreto confinado
c       euc=deformacion al final de la rama descendente del concre-
c           to  no confinado
c       eum=deformacion maxima del concreto confinado
c       cc=resultante de las fuerzas de compresion del concreto
c       mc=momento resistente debido a las fuerzas de compresion
c          del concreto
c       eccu=deformaciÛn ultima del concreto confinado
c
c******************************************************************************
c
	real*8 ec,hc,x,h,fc,eo,rv,k
	real*8 euc,eum,eccu,b
	real*8 cc1,mc1
c**********************************************************************
c       LOCALES
c**********************************************************************
	real*8 ecc,Yn,c1a,m1a,c1aa,c2aa,c1b,m1b,c1bb,c2bb,q
	real*8 Ccnc,Mcnc,ecco
	real*8 comp1,comp2,comp3,brazo1,brazo2,brazo3,yu
c
c*******************************************************************
c      Calculo de la fuerza resultante de compresion del concreto de
c      la parte NO CONFINADA (HOGNESTAD)
c*******************************************************************
        ecc=ec
        if ((ec.gt.euc).and.(ecc.gt.euc)) then
	      Yn=x*euc/ec
            call f_comp_concret(euc,b,Yn,h,fc,eo,euc,
     s                             c1a,m1a,c1aa,c2aa)
            call f_comp_concret(euc,hc,Yn,h,fc,eo,euc,

     s                             c1b,m1b,c1bb,c2bb)
            q=(ec-euc)/euc
            Ccnc=c1a-c1b
            Mcnc=(m1a-(c1aa+c2aa)*(q*Yn+rv))-(m1b-(c1bb+c2bb)*(q*Yn+rv))
	  end if
c       Calculo de la fuerza resultante de compresion del concreto de
c	  la parte  CONFINADA (Kent y Park)
c

	  ecco=k*0.002
	  if ((ecc.ge.eccu).and.(ecc.le.eum)) then
                  comp1=2.*k*fc*hc*(x)*(ecco/ecc)/3.
	            brazo1=5.*ecco*(x-rv)/(8.*ecc)
	            comp2=hc*(x)*(1.2*k*fc)*(eccu-ecco)/(2*eccu)
	            brazo2=(x)*((ecco/eccu)+((1.-ecco/eccu)*(2.2*k*fc))
     1                      /(3.*1.2*k*fc))
	            yu=(x)*eccu/ecc
	            comp3=0.2*k*fc*hc*(x-yu)
              	brazo3=0.5*(x+yu)
	  endif
c***********************************************************************
c
c       Calculo del momento resistente debido a las fuerzas de compresion 
c           del
c	  concreto
c***********************************************************************
c
	  mc1=comp1*(h/2.-x-rv+brazo1)+comp2*(h/2.-x-rv+brazo2)+comp3*
     1             (h/2.-x-rv+brazo3)+Mcnc
	  cc1=comp1+comp2+comp3+Ccnc
c
	return
	end
c
c
c*********************************************************************************************
c        SUBRUTINA INTERPOLACION
c        INTERPOLACION DE UNA TABLA DE VALORES QUE DEFINEN UN 
c        DIAGRAMA DE INTERACCION
c*********************************************************************************************
c
        subroutine interpolacion(p,max_tabla,tabla,n_tabla,
     s                                           var)
c
       IMPLICIT NONE
c*********************************************************************************************
c        DECLARACIONES
c*********************************************************************************************
c
        integer max_tabla,i,n_tabla
        real*8  p,var,tabla(max_tabla,2)
c
c*********************************************************************************************
c        VERIFICACION DE LOS VALORES DE LAS TABLAS
c        (los valores de las fuerzas axiales deben ir de mayor a menor)
c*********************************************************************************************
c
        do i=1,n_tabla-1
                if(tabla(i+1,1).ge.tabla(i,1))then
                        stop "error en los datos de los diagramas"
                end if
        end do
c
c*********************************************************************************************
c
        if((p.gt.tabla(1,1)).or.(p.lt.tabla(n_tabla,1)))then
                var=0.
        else
                i=1
                do while((p.le.tabla(i,1)).and.(i.le.n_tabla))
                        i=i+1
                        if(p.ge.tabla(i,1))then
                         var=tabla(i,2)+
     1                   (tabla(i,2)-tabla(i-1,2))*(p-tabla(i,1))/
     2                   (tabla(i,1)-tabla(i-1,1))
                        end if
                end do
      end if
c
      return
      end  
c
c*********************************************************************************************
c
c******************************************************************************
c       SUBRUTINA INTERPOLACION
c       INTERPOLACION DE UNA TABLA DE VALORES QUE DEFINEN UN
c       DIAGRAMA DE INTERACCION
c
c******************************************************************************
c
	subroutine interpolacion_x(p,max_tabla,tabla,n_tabla,
     s                          var)
c
c******************************************************************************
c       DECLARACIONES
c******************************************************************************
c
	integer max_tabla,i,n_tabla
	real*8 tabla(max_tabla,3),var,p
c
c*****************************************************************************
c       VERIFICACION DE LA TABLA POR SI ACASO
c       (los valores de las fuerzas axiales deben ir de menor a mayor)
******************************************************************************
c
	do i=1,n_tabla-1
		if(tabla(i+1,1).ge.tabla(i,1))then
			stop "error en los datos de los diagramas"
		end if
	end do
c
c******************************************************************
c
	if((p.gt.tabla(1,1)).or.(p.lt.tabla(n_tabla,1)))then
		var=0.
	else
		i=1
		do while((p.le.tabla(i,1)).and.(i.le.n_tabla))
			i=i+1
			if(p.ge.tabla(i,1))then
				var=tabla(i,3)+
     1                    (tabla(i,3)-tabla(i-1,3))*(p-tabla(i,1))/
     2                    (tabla(i,1)-tabla(i-1,1))
			end if
		end do
	end if
c
	return
	end
c
c******************************************************************
c
c       SUBRUTINA CAL_COEFIC PARA CUALQUIER MODELO
c       CALCULO DE LOS COEFICIENTES c,q,gcr,my DEL MODELO DE PORTI-
c       COS ELASTO-PLASTICOS  DEGRADABLES  A PARTIR DEL MOMENTO
c       PLASTICO, ULTIMO, DE AGRIETAMIENTO Y DE LA ROTULA PLASTICA
c       ULTIMA  MEDIANTE    LA    RESOLUCION    DE   DOS    SISTEMAS  DE
c       ECUACIONES NO  LINEALES  A PARTIR DE LAS FUNCIONES DE DANO
c       FLUENCIA
c
c******************************************************************
c
        subroutine cal_coefic(mcr,mp,mu,fpu,ei,l,
     s                       gcr,q,m,c,dp,du)   
c
       IMPLICIT NONE
c******************************************************************
c               DECLARACIONES LOCALES
c******************************************************************
c
c       mcr= Valor critico o de agrietamiento
c       mp= Valor plastico
c       mu= Valor ultimo
c       fpu=rotula plastica ultima
c       ei=modulo de elasticidad*inercia
c       l=longitud del elemento
c       gcr=coeficiente del criterio de plasticidad
c       q=coeficiente del criterio de plasticidad
c       my=coeficiente del criterio de plasticidad
c       c=coeficiente del criterio de plasticidad
c       f0=coeficiente de flexibilidad a la flexion
c       fp=rotacion plastica
c       m=momento
c       d=variable dano
c       dp=valor de la variable dano para el momento plastico
c       du=valor de la variable dano para el momento ultimo
c       g=funcion de dano
c       f=funcion de fluencia
c       ff=derivada parcial del momento con respecto al dano
c
c******************************************************************
c
        real*8 mcr,mp,mu,fpu
        real*8 ei,l,f0
        real*8 gcr,q,c,m
        real*8 dmin,dmax,dp,du
        real*8 fmin,fmax,fu
        real*8 gmin,gmax,gp

c
c******************************************************************
c       Calculo del coeficiente de flexibilidad a flexion
c
             f0=l/(3.*ei)
c
c       Calculo del coeficiente "gcr" del criterio de plasticidad
c       para m=mcr, d=0 y g=0
c
              gcr=(mcr*mcr)*f0/(2.)
c
c       Calculo del coeficiente "q" del criterio de plasticidad 
c       y de du para m=mu, g=0 y ff=0
c       Resolucion del sistema de ecuaciones no lineales para el
c       calculo de "du" de la seccion y del coeficiente "q"
              dmin=0.0005
              dmax=0.9995
              call func_f(f0,gcr,mu,dmin,fmin,q)
              call func_f(f0,gcr,mu,dmax,fmax,q)
c
        do while((dmax-dmin).gt.0.00001)
              du=(dmin+dmax)/2.
              call func_f(f0,gcr,mu,du,fu,q)
c
          if(fu.gt.0) then
              dmin=du
              fmin=fu
          else
              dmax=du
              fmax=fu
          end if
c
        end do
c
c       Calculo del coeficiente "my" del criterio de plasticidad y
c       de dp para g=0, f=0, m=mp, fp=0
c       Calculo del coeficiente "c" del criterio de plasticidad 
c       para m=mu, d=du, f=0 y fp=fpu
c       Resolucion de la ecuacion para el calculo de dp de la sec-
c       cion y los coeficientes c y my
              dmin=0         
              dmax=du
              call func_g(f0,gcr,mp,dmin,q,gmin)
              call func_g(f0,gcr,mp,dmax,q,gmax)
c
        do while((dmax-dmin).gt.0.00001)
                dp=(dmin+dmax)/2.
                call func_g(f0,gcr,mp,dp,q,gp)
c
          if(gp.gt.0) then
                dmin=dp
                gmin=gp
          else
                dmax=dp
                gmax=gp
          end if
c
        end do
c
c
              m = mp/(1.-dp)
c
c
              c = (mu/(1-du) - m)/fpu
c
c
        return
        end
c
c******************************************************************
c
c       SUBRUTINA func_f PARA CUALQUIER MODELO SE EVALUA LA FUNCION
c       DE DANO Y  LA  DERIVADA  DEL  MOMENTO CON RESPECTO AL DA•O
c       PARA RESOLVER EL SISTEMA DE ECUACIONES NO LINEALES PARA LA
c       DETERMINACION DEL COEFICIENTE q Y DE du
c
c******************************************************************
c
        subroutine func_f(f0,gcr,m,d,ff,q)
c
       IMPLICIT NONE
c******************************************************************
c
        real*8 f0,gcr,m,d,ff,q
c
c******************************************************************
c       
              q=((m*m)*f0/(2.)-gcr*(1.-d)*(1.-d))/((1.-d)*dlog(1.-d))
c
              ff=-2.*gcr*(1.-d)-q*(dlog(1.-d)+1.)
c
        return
        end
c
c******************************************************************
c
c       SUBRUTINA func_g PARA CUALQUIER MODELO SE EVALUA LA FUNCION
c       DE DANO PARA RESOLVER EL SISTEMA DE ECUACIONES NO LINEALES
c       PARA LA DETERMINACION DEL COEFICIENTE my  Y  DE  dp
c
c******************************************************************
c
        subroutine func_g(f0,gcr,m,d,q,g)
c
       IMPLICIT NONE
c******************************************************************
c
        real*8 f0,gcr,m,d,g,q
c
c******************************************************************
c
              g=(m*m)*f0/(2.)-gcr*(1.-d)*(1.-d)-q*(1.-d)*dlog(1.-d)
c
        return
	  end
c
c******************************************************************
c
c      CALCULO DEL DA•O POR CORROSION
c      ESCRITA POR: Scarlet Montilla, Carlos, Julio Fl¢rez, Ricardo Pic¢n
c      A•O: 2021
c
c******************************************************************
c
       subroutine cal_cor(
     e            fp,d,dtiempo,tiempo_total,t_ini,phi,io,
     e            iod_I,dI,Rel,corr_cero,DhDq,
     s            corr,fccorr)
c
c
       IMPLICIT NONE
c       t_ini: Tiempo inicial de corrosi¢n
c       phi: Promedio de los di†metros de las barras
c       io: Densidad inicial de corriente de corrosi¢n
c       iod_I: Densidad de corriente de identificaci¢n
c       dI: Da§o de identificaci¢n
c       Rel: Par†metro Probabil°stico
c       cor: Da§o en la r¢tula por corrosi¢n
c       fcc: Fuerza termodin†mica en la r¢tula por corrosi¢n
c       DhDq: ( @h/@c / @q/@c )
c       fp: Deformacion plastica de la r¢tula
c       d: da§o de la r¢tula
c       dtime: Intervalo de tiempo
c       tiempo_total: Tiempo total de an†lisis
c
       real*8 t_ini,phi,io,iod_I,dI,Rel,corr,fccorr,tiempo_total,
     1        dtiempo,corr_cero,fp,d,DhDq,aux1,aux2
c
       real*8 o_qui,ttr,delta_tiempo,o_mec
c
       t_ini=t_ini/3.154e7
	 ttr = tiempo_total/3.154e7
	 delta_tiempo = dtiempo/3.154e7
c
       if (ttr.le.t_ini) then
            corr = 0.
            fccorr = 0.
       else
           if (ttr.le.(t_ini+.57)) then
             o_qui = 0.06*io/phi
           else
             o_qui = (0.036*io*(ttr-t_ini)**(-0.29))/phi
           end if
           aux1 = (0.0116*Rel*(iod_I-io))/(phi*(dlog(1.-dI))**2)
           aux2 = (dlog(1.-d))**2-fp**2*(1.-d)*DhDq
           o_mec = aux1 * aux2
           if (o_mec.le.0.) then 
		     o_mec= 0.0
		   end if
		   write(*,*)DhDq,o_qui,o_mec,aux1,aux2,iod_I,io,Rel
	
           corr = corr_cero +(o_qui+o_mec)*delta_tiempo
           fccorr = 0.
       end if
c
c
      return
      end
c
c
c*******************************************************************
c
c        Propiedades de las Barras de acero por corrosi¢n
c
c*******************************************************************
c
c
       subroutine Prop_acero_i(
     e            prop,corr_cero)
c
c
       IMPLICIT NONE
      Real(kind=8):: corr_cero,prop(61),Es,kcor,aux
c
c*******************************************************************************************************************
c       MODIFICACI‡N DE LAS PROPIEDADES DEL ACERO
c        POR CORROSI‡N
c******************************************************************
         Es=prop(11)/prop(13)
c
c        Extremo i
c
            prop(11)=(1.-0.1*corr_cero)*prop(11)
            prop(12)=(1.-0.15*corr_cero)*prop(12)
            prop(13)=prop(11)/Es
            prop(14)=(1.-0.2*corr_cero)*prop(14)
            prop(15)=(1.-0.25*corr_cero)*prop(15)
c
c********************************************************************
c      REDUCCI‡N DEL AREA EFECTIVA DE LAS BARRAS DE ACERO
c                         POR CORROSION
c                     Acero en el extremo "i"
c                 Prop 17, 19, 21, 23, 25, 27, 29
c********************************************************************
c
            kcor=-4.*(corr_cero)**2.*
     1         (asin(sqrt(-(corr_cero)**2.+1.)))+
     2         2.*corr_cero*(sqrt(-(corr_cero)**2.+1.))
c
        if (corr_cero.lt.(sqrt(2.)/2.)) then
           aux=kcor+3.1415926535898-(asin(2.*corr_cero*
     1         (sqrt(-(corr_cero)**2.+1.))))
c
           prop(17)=prop(17)*aux/3.1415926535898
           prop(19)=prop(19)*aux/3.1415926535898
           prop(21)=prop(21)*aux/3.1415926535898
           prop(23)=prop(23)*aux/3.1415926535898
           prop(25)=prop(25)*aux/3.1415926535898
           prop(27)=prop(27)*aux/3.1415926535898
           prop(29)=prop(29)*aux/3.1415926535898
        else
           aux=kcor+(asin(2.*corr_cero*
     1         (sqrt(-(corr_cero)**2.+1.))))
c
           prop(17)=prop(17)*aux/3.1415926535898
           prop(19)=prop(19)*aux/3.1415926535898
           prop(21)=prop(21)*aux/3.1415926535898
           prop(23)=prop(23)*aux/3.1415926535898
           prop(25)=prop(25)*aux/3.1415926535898
           prop(27)=prop(27)*aux/3.1415926535898
           prop(29)=prop(29)*aux/3.1415926535898
c
        end if
      return
      end
c
c********************************************
c          Extremo j
c********************************************
c
       subroutine Prop_acero_j(
     e            prop,corr_cero)
c
       IMPLICIT NONE
c
      Real(kind=8):: corr_cero,prop(61),Es,kcor,aux
c
c*******************************************************************************************************************
c       MODIFICACI‡N DE LAS PROPIEDADES DEL ACERO
c        POR CORROSI‡N
c******************************************************************
         Es=prop(11)/prop(13)
c
         prop(11)=(1.-0.41*corr_cero)*prop(11)
         prop(12)=(1.-0.45*corr_cero)*prop(12)
         prop(13)=prop(11)/Es
         prop(14)=(1.-0.8*corr_cero)*prop(14)
         prop(15)=(1.-1.1*corr_cero)*prop(15)
c
c
c********************************************************************
c      REDUCCI‡N DEL AREA EFECTIVA DE LAS BARRAS DE ACERO
c                         POR CORROSION
c                     Acero en el extremo "j"
c                Prop 41, 43, 45, 47, 49, 51, 53
c********************************************************************
c
            kcor=-4.*(corr_cero)**2.*
     1         (asin(sqrt(-(corr_cero)**2.+1.)))+
     2         2.*corr_cero*(sqrt(-(corr_cero)**2.+1.))
c
        if (corr_cero.lt.(sqrt(2.)/2.)) then
           aux=kcor+3.1415926535898-(asin(2.*corr_cero*
     1         (sqrt(-(corr_cero)**2.+1.))))
c
           prop(41)=prop(41)*aux/3.1415926535898
           prop(43)=prop(43)*aux/3.1415926535898
           prop(45)=prop(45)*aux/3.1415926535898
           prop(47)=prop(47)*aux/3.1415926535898
           prop(49)=prop(49)*aux/3.1415926535898
           prop(51)=prop(51)*aux/3.1415926535898
           prop(53)=prop(53)*aux/3.1415926535898
        else
           aux=kcor+(asin(2.*corr_cero*
     1         (sqrt(-(corr_cero)**2.+1.))))
c
           prop(41)=prop(41)*aux/3.1415926535898
           prop(43)=prop(43)*aux/3.1415926535898
           prop(45)=prop(45)*aux/3.1415926535898
           prop(47)=prop(47)*aux/3.1415926535898
           prop(49)=prop(49)*aux/3.1415926535898
           prop(51)=prop(51)*aux/3.1415926535898
           prop(53)=prop(53)*aux/3.1415926535898
c
        end if
c
c
      return
      end
c