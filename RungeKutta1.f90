!--------------------------
!RungeKuttaMethod
!
!dy/dx=F(x,y)
!
!20150820
!--------------------------

program main
implicit none

double precision,parameter::xini=0.0d0		!x初期値
double precision,parameter::yini=0.50d0		!y初期値
double precision,parameter::xfin=50.0d0		!xどこまで計算するか
integer,parameter::xNbin=100				!刻み数

integer::i
double precision::x,y
double precision::dx
double precision::k1,k2,k3,k4
open(10,file='runge_kutta2.d')

dx=(xfin-xini)/dble(xNbin)		!刻む幅を計算
x=xini 							!初期値を代入
y=yini							!初期値を代入

!ルンゲクッタ法
do i=1,xNbin
	k1=dx*func(x,y)
	k2=dx*func(x+0.50d0*dx,y+0.50d0*k1)
	k3=dx*func(x+0.50d0*dx,y+0.50d0*k2)
	k4=dx*func(x+dx,y+k3)
	y=y+(k1+2d0*k2+2d0*k3+k4)/6d0
    x=xini+dx*dble(i)
    write(10,'(2e26.16)') x,y
end do

close(10)

!-----------------------------------------

stop
contains

!解きたい微分方程式

function func(X,Y) result(dydx)
	double precision,intent(in)::X,Y
	double precision dydx

	dydx=y*(1-y)  		!この式を解く

end function func

end program main