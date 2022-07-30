    subroutine eigenmatrix(R,L,rho,u,v,p)
    
    integer dim
    parameter(dim = 4)
    real R(dim,dim),L(dim,dim),rho,u,v,p,c,H,VV,gamma,gamma1,E
    
    gamma = 1.4d0
    
    gamma1 = gamma - 1
    
    E = p/gamma1 + 0.5d0*rho*(u**2 + v**2)
    
    c = (abs(gamma*p/rho))**0.5d0
    H = (E + p)/rho
    VV = u**2 + v**2
            
    ! 组装右特征向量
    R(1,1) = 1
    R(2,1) = u - c
    R(3,1) = v
    R(4,1) = H - u*c
            
    R(1,2) = 1
    R(2,2) = u
    R(3,2) = v
    R(4,2) = 0.5d0*VV
            
    R(1,3) = 0
    R(2,3) = 0
    R(3,3) = 1
    R(4,3) = v
            
    R(1,4) = 1
    R(2,4) = u + c
    R(3,4) = v
    R(4,4) = H + u*c
            
    ! 组装左特征向量
    L(1,1) = H + (c/gamma1)*(u - c)
    L(1,2) = -(u + c/gamma1)
    L(1,3) = -v
    L(1,4) = 1
            
    L(2,1) = -2*H + 4d0*c**2/gamma1
    L(2,2) = 2*u
    L(2,3) = 2*v
    L(2,4) = -2
            
    L(3,1) = -2*v*c**2/gamma1
    L(3,2) = 0
    L(3,3) = 2*c**2/gamma1
    L(3,4) = 0
            
    L(4,1) = H - (c/gamma1)*(u + c)
    L(4,2) = -u + c/gamma1
    L(4,3) = -v
    L(4,4) = 1
            
    L = (gamma1/(2*c**2))*L
    
    ! 不做特征分解
    !L = 0
    !R = 0
    !L(1,1) = 1
    !L(2,2) = 1
    !L(3,3) = 1
    !L(4,4) = 1
    !R(1,1) = 1
    !R(2,2) = 1
    !R(3,3) = 1
    !R(4,4) = 1
    
    end subroutine eigenmatrix