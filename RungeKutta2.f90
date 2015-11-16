!-------------------
!runge_kutta_method
!
!funcdy:dy/dx=F1(x,y,z)
!dz/dx=F2(x,y,z)
!
!20150820
!-------------------

program main
implicit none

double precision,parameter::xini=0.0d0		!x初期値
double precision,parameter::yini=0.0d0		!y初期値
double precision,parameter::zini=1.0d0		!z初期値

double precision,parameter::xfin=10.0d0		!xをどこまで計算するか
!double precision,parameter::zfin=50.0d0
integer,parameter::xNbin=10000				!刻み数

integer::i
double precision::x,y,z
double precision::dx
double precision::k1,k2,k3,k4,q1,q2,q3,q4
open(10,file='runge_kutta2.d')

dx=(xfin-xini)/dble(xNbin)		!刻む幅を計算
x=xini 							!初期値を代入
y=yini							!初期値を代入
z=zini

!ルンゲクッタ法の計算をする
do i=1,xNbin
	k1=dx*funcdy(x,y,z)
	q1=dx*funcdz(x,y,z)

	k2=dx*funcdy(x+0.50d0*dx,y+0.50d0*k1,z+0.50d0*q1)
	q2=dx*funcdz(x+0.50d0*dx,y+0.50d0*k1,z+0.50d0*q1)
	
	k3=dx*funcdy(x+0.50d0*dx,y+0.50d0*k2,z+0.50d0*q2)
	q3=dx*funcdz(x+0.50d0*dx,y+0.50d0*k2,z+0.50d0*q2)
	
	k4=dx*funcdy(x+dx,y+k3,z+q3)
	q4=dx*funcdz(x+dx,y+k3,z+q3)

	y=y+(k1+2d0*k2+2d0*k3+k4)/6d0
	z=z+(q1+2d0*q2+2d0*q3+q4)/6d0
    x=xini+dx*dble(i)
    write(10,'(2e26.16)') x,y,z
end do

close(10)

!-----------------------------------------

stop
contains


!dy/dx
function funcdy(X,Y,Z) result(dydx)		
	double precision,intent(in)::X,Y,Z
	double precision dydx

	dydx=Z  		

end function funcdy

!dz/dx
function funcdz(X,Y,Z) result(dzdx)
	double precision,intent(in)::X,Y,Z
	double precision dzdx

	dzdx=-Y  		

end function funcdz

end program main