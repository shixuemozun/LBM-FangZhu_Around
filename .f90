module para
integer,parameter::lx=600,ly=300
real*8::U11 = 0.1
real*8::ex(0:8)=[0.d0,1.d0,0.d0,-1.d0,0.d0,1.d0,-1.d0,-1.d0,1.d0],&
        ey(0:8)=[0.d0,0.d0,1.d0,0.d0,-1.d0,1.d0,1.d0,-1.d0,-1.d0]
real*8::w(0:8)=[4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36]
real*8::rho(0:lx,0:ly),u(0:lx,0:ly),v(0:lx,0:ly),ff(0:lx,0:ly,0:8),f(0:lx,0:ly,0:8),feq(0:lx,0:ly,0:8)
real*8::P(0:lx,0:ly)
integer::i,j,k,ip,jp,n,x,y
real*8::c,cs,Re,dx,dy,lx_,ly_,dt,rho0,T0,P0,tau_f,niu,error
real*8::eu,uv,rho1,u0,u1
real*8::flag(0:lx,0:ly)
integer::t1 = 0

end module para




program main
use para
implicit none
call system('mkdir output')


dx = 1.0
dy = 1.0
lx_ = dx*dble(lx)
ly_ = dy*dble(ly)
dt = dx
c = dx/dt
rho0 = 1.0
Re = 500
niu = (U11*lx_)/Re
tau_f = 3.0*niu+0.5


do i=0,lx
do j=0,ly
flag(i,j) = 0

enddo
enddo

!设定区域标号


do i=99,129
do j=124,184
flag(i,j) = 1
enddo
enddo
!方柱信号标记



do i=0,lx
do j=0,ly
u(i,j) = 0
v(i,j) = 0
rho(i,j) = 1.0
u(0,j) = -1*(DBLE(j)/(4*ly_))*(DBLE(j)/(4*ly_))+0.25*(DBLE(j)/(4*ly_))

enddo
enddo


do i=0,lx
do j=0,ly
P(i,j) = 0

enddo
enddo


do i=0,lx
do j=0,ly
do k=0,8
eu = ex(k)*u(i,j)+ey(k)*v(i,j)
uv = u(i,j)**2+v(i,j)**2
feq(i,j,k) = w(k)*rho(i,j)*(1.0+3.0*eu+4.5*eu**2-1.5*uv)
enddo
enddo
enddo

do i=0,lx
do j=0,ly
do k=0,8
f(i,j,k) = feq(i,j,k)
enddo
enddo
enddo

t1 = 0
do while(t1<5000)
t1 = t1+1
write(*,'(I10)')t1

do i=1,lx-1
do j=1,ly-1
do k = 0,8
ip = i-ex(k)
jp = j-ey(k)

ff(i,j,k) = f(ip,jp,k)+(feq(ip,jp,k)-f(ip,jp,k))/1.2
enddo
enddo
enddo





do i=1,lx-1
do j=1,ly-1
do k=0,8
f(i,j,k) = ff(i,j,k)
enddo

rho1 = 0
do k=0,8
rho1 = rho1+f(i,j,k)
enddo
rho(i,j) = rho1

u0 = 0
do k=0,8
u0 = u0 +ex(k)*f(i,j,k)
enddo
u(i,j) = u0/rho(i,j)

u1 = 0
do k=0,8
u1= u1+ ey(k)*f(i,j,k)
enddo
v(i,j) = u1/rho(i,j)
enddo
enddo

do i=99,129
rho(i,124) = rho(i,123)
u(i,124) = 0
v(i,124) = 0
do k=0,8
f(i,124,k) = feq(i,124,k)+f(i,123,k)-feq(i,123,k)
enddo
rho(i,184) = rho(i,185)
u(i,184) = 0
v(i,184) = 0
do k=0,8
f(i,184,k) = feq(i,184,k)+f(i,185,k)-feq(i,185,k)
enddo
enddo

!柱体上下


do j=125,183
rho(99,j) = rho(98,j)
u(99,j) = 0
v(99,j) = 0
do k=0,8
f(99,j,k) = feq(99,j,k)+f(98,j,k)-feq(98,j,k)
enddo

rho(129,j) = rho(130,j)
u(129,j) = 0
v(129,j) = 0
do k=0,8
f(129,j,k) = feq(129,j,k)+f(130,j,k)-feq(130,j,k)
enddo
enddo
!柱体左右


do j=1,ly-1
do k=0,8
rho(lx,j) = rho(lx-1,j)
rho(0,j) = rho(1,j)

u(0,j) = -1*(DBLE(j)/(4*ly_))*(DBLE(j)/(4*ly_))+0.25*(DBLE(j)/(4*ly_))
v(0,j) = 0
u(lx,j) = u(lx-1,j)
v(lx,j) = v(lx-1,j)
f(0,j,k) = feq(0,j,k)+f(1,j,k)-feq(1,j,k)
f(lx,j,k) = feq(lx,j,k)+f(lx-1,j,k)-feq(lx-1,j,k)

enddo
enddo


do i=0,lx
do k=0,8

rho(i,0) = rho(i,1)
f(i,0,k) = feq(i,0,k)+(f(i,1,k)-feq(i,1,k))
u(i,0) = 0
v(i,0) = 0
u(i,ly) = 0
v(i,ly) = 0
rho(i,ly) = rho(i,ly-1)
f(i,ly,k) = feq(i,ly,k)+f(i,ly-1,k)-feq(i,ly-1,k)
enddo
enddo



do i=0,lx
do j=0,ly
do k=0,8
eu = ex(k)*u(i,j)+ey(k)*v(i,j)
uv = u(i,j)**2+v(i,j)**2
feq(i,j,k) = w(k)*rho(i,j)*(1.0+3.0*eu+4.5*eu**2-1.5*uv)
enddo
enddo
enddo

if(mod(t1,500)==0) then
open(100,file='./output/Flow.dat')

write(100,*) 'variables = x,y,u,v'
write(100,*) 'zone j=',lx+1,',k=',ly+1,',f=point'

do y=0,ly
do x =0,lx

write(100,'(10D20.10)')dble(x)/lx_,dble(y)/(4*ly_),u(x,y),v(x,y)
enddo
enddo
close(100)
endif






enddo !while的结束语句





end !主程序的结束句