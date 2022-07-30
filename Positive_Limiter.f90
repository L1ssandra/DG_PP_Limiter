    subroutine Positive_Limiter(Pol,PolValue,hx,hy,dt,Nx,Ny,dim,dimpol,diagpol)
    
    integer i,j,Nx,Ny,dim,dimpol,quad,i1,j1,d1,check,k1
    parameter(quad = 4)
    real Pol(Nx,Ny,dimpol,dim),PolValue(quad,quad,dimpol,dim,2),PolDx(quad,quad,dimpol,dim),PolDy(quad,quad,dimpol,dim)
    real Poledge(4,quad,dimpol,dim),Q(Nx,Ny,dim),diagpol(Nx,Ny,dimpol,dim),DPol(Nx,Ny,dimpol,dim),alphax,alphay,alphatest
    real PolR(Nx,Ny,dimpol,dim),PolL(Nx,Ny,dimpol,dim),PolU(Nx,Ny,dimpol,dim),PolD(Nx,Ny,dimpol,dim)
    real u,v,p,rho,gamma
    real IV(dim),IR(dim),IL(dim),IU(dim),ID(dim),lambdai(quad),lambdaj(quad),weight(quad),fu(dim,2)
    real uint(dim),uext(dim),hLF(dim),fuint(dim),fuext(dim),hx1,hy1,minrho,minrho0
    real Poltilde(Nx,Ny,dimpol,dim),ta,thetamin,sa(dim)
    
    check = 0
    minrho0 = 1e-13
    
    gamma = 1.4d0
    
    ! 先限制rho
    do i = 1,Nx
        do j = 1,Ny
            
            minrho = Pol(i,j,1,1)
            !if (Pol(i,j,1,1) < minrho0) then
            !    print *,Pol(i,j,1,1)
            !    Pol(i,j,1,1) = minrho0
            !    Pol(i,j,2:dimpol,1) = 0
            !end if
            
            do i1 = 1,quad
                do j1 = 1,quad
                    do k1 = 1,2
                        uint = 0
                        do d1 = 1,dimpol
                            uint = uint + Pol(i,j,d1,:)*PolValue(i1,j1,d1,:,k1)
                        end do
                        if (uint(1) < minrho) then
                            minrho = uint(1)
                        end if
                    end do
                end do
            end do
            
            theta = min( abs((minrho0 - Pol(i,j,1,1))/(minrho - Pol(i,j,1,1))) , 1d0 )
            !if (theta < 1) then
            !    print *,"改变了密度，theta = ",theta
            !end if
            Pol(i,j,2:dimpol,1) = theta*Pol(i,j,2:dimpol,1)
        end do
    end do
    
    minrho0 = 1e-13
    ! 再限制p
    do i = 1,Nx
        do j = 1,Ny
            
            check = 0
            thetamin = 1d0
            
            do i1 = 1,quad
                do j1 = 1,quad
                    do k1 = 1,2
                        uint = 0
                        do d1 = 1,dimpol
                            uint = uint + Pol(i,j,d1,:)*PolValue(i1,j1,d1,:,k1)
                        end do
                        if (uint(1) < 0) then
                            !print *,uint(1)
                        end if
                    
                        if ( uint(4) - 0.5d0*(uint(2)**2 + uint(3)**2)/uint(1) < 0) then
                            !print *,"积分点",i1,j1,"的值为",uint(4) - 0.5d0*(uint(2)**2 + uint(3)**2)/uint(1)
                            !print *,"单元平均值为",Pol(i,j,1,4) - 0.5d0*(Pol(i,j,1,2)**2 + Pol(i,j,1,3)**2)/Pol(i,j,1,1)
                            !print *,i1,j1,k1
                            !check = 1
                        end if
                        if ( uint(4) - 0.5d0*(uint(2)**2 + uint(3)**2)/uint(1) >= minrho0/(gamma - 1) ) then
                            ta = 1
                        else
                            !print *,Pol(i,j,1,:),uint
                            call Calculate_ta(ta,Pol(i,j,1,:),uint)
                        end if
                        if (ta < 1) then
                            !print *,"ta = ",ta
                        end if
                    
                        sa = ta*uint + (1 - ta)*Pol(i,j,1,:)
                    
                        call Calculate_theta(theta,sa,Pol(i,j,1,:),uint)
                    
                        if (theta < thetamin) then
                            thetamin = theta
                        end if
                    
                    end do
                end do
            end do
            
            if (thetamin < 1) then
                !print *,"改变了压强，theta = ",thetamin
                !check = 1
                !thetamin = 0.99d0*thetamin
            end if
            
            Pol(i,j,2:dimpol,:) = thetamin*Pol(i,j,2:dimpol,:)
            
            ! 检测一下
            if (check == 1) then
                print *,"开始检测↓"
                do i1 = 1,quad
                    do j1 = 1,quad
                        do k1 = 1,2
                            uint = 0
                            do d1 = 1,dimpol
                                uint = uint + Pol(i,j,d1,:)*PolValue(i1,j1,d1,:,k1)
                            end do
                            !if (uint(4) - 0.5d0*(uint(2)**2 + uint(3)**2)/uint(1) < minrho0) then
                            !    print *,uint(4) - 0.5d0*(uint(2)**2 + uint(3)**2)/uint(1)
                            !end if
                        end do
                    end do
                end do
                print *,"检测结束↑"
            end if
            
        end do
    end do
    
    end subroutine Positive_Limiter
    
    
    
    subroutine Calculate_theta(theta,s,w,q)
    
    real s(4),w(4),q(4),theta,sw(4),qw(4)
    
    sw = s - w
    qw = q - w
    
    
    theta = (sw(1)**2 + sw(2)**2 + sw(3)**2 + sw(4)**2)**0.5d0/(qw(1)**2 + qw(2)**2 + qw(3)**2 + qw(4)**2)**0.5d0
    
    end subroutine Calculate_theta
    
    
    subroutine Calculate_ta(ta,w,q)
    
    real ta,w(4),q(4),a(4),b(4),gamma,epsilon,epsilon1,t1,t2,tc,uc(4),p
    
    a = q - w
    b = w
    
    epsilon = 1e-6
    gamma = 5d0/3d0
    
    t1 = 0
    t2 = 1
    
    do while (t2 - t1 > epsilon)
        tc = 0.5d0*(t1 + t2)
        uc = a*tc + b
        call pressure(uc,p) 
        if (p > epsilon) then
            t1 = tc
        else
            t2 = tc
        end if
    end do
    
    ta = t1
    
    end subroutine Calculate_ta
    
    
    
    subroutine pressure(u,p)
    
    real u(4),p,gamma1
    
    gamma1 = 1.4d0
    
    p = gamma1*(u(4) - 0.5d0*(u(2)**2 + u(3)**2)/u(1))
    
    end subroutine pressure
    
    