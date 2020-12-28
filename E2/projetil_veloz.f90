!Programa que calcula a velocidade de um projétil durante sua trajetória

	program projetil_veloz
	IMPLICIT NONE

	!Variáveis
	real*8 :: vo !velocidade inicial
	real*8 :: vx !velocidade em x
	real*8 :: vy !velocidade em y
	real*8 :: y !posição y
	real*8 :: x !posição x
	real*8 :: teta !angulo da velocidade com o eixo x
	real*8 :: m !massa
	real*8 :: g !gravidade
	real*8 :: r !densidade
	real*8 :: ro !densidade inicial
	real*8 :: A !área
	real*8 :: CD_menor !drag coefficient para velocidades abaixo da velocidade do som
	real*8 :: CD_maior !drag coefficient para velocidades acima da velocidade do som
	real*8 :: dt !Variação de tempo
	real*8 :: aux_x !auxiliar para x
	real*8 :: aux_y !auxiliar para y
	real*8 :: aux
	real*8 :: v_som !velocidade do som
	real*8 :: beta !taxa de variação da temperatura com a altitude
	real*8 :: alpha !expoente para o ar
	real*8 :: To !temperatura a nível do mar
	real*8 :: cal !calibre
	real*8 :: v !velocidade
	real*8 :: t !tempo
	real, parameter :: pi = 3.1415927

	!Inicalizando variáveis:
	vo = 373.23d0 !m/s, dado
	m = 42.d0 !kg, dado
	cal = 149.1d0 !mm, dado
	g = 9.8d0 !m/s², dado
	CD_menor = 0.295d0 !dado
	CD_maior = 0.5d0 !dado
	v_som = 343.d0 !dado
	alpha = 4.256d0 !dado
	beta = 6.43d-3 !dado
	To = 300 !K, dado
	ro = 0.d0 !kg/m³, dado
	dt = 0.01d0 !dado
	teta = 0.74267252397537231d0 !obtido com o programa alcance_res
	vx = vo*(cos(teta))
	vy = vo*(sin(teta))
	x = 0.d0
	y = 0.d0
	t = 0.d0
	v = vo

	!Calculando a área:
	A = pi*(((cal*0.001d0)/2.d0)**2.d0) !pi*raio²

	!Abrindo arquivo para impressão de sesultados
	!Aonde escreveremos os alcances
	open(36,file = "projetil_veloz_2.dat")
	write(36,*)"teta: ", teta," velocidade: ", vo
	!aonde escreveremos dados para gráficos
	open(37,file = "projetil_veloz_2_graph.dat")


	!Calculamos o movimento de lançamento assumindo que o objeto começa na mesma altura que termina, ou seja, zero.

	DO WHILE(y >= 0)

		!Utilizamos o método de Euler para calcular a altura e alcance do objeto a cada pedaço de tempo, e imprimimos a trajetória.

		aux_x = x
		aux_y = y

		x = aux_x + vx*dt
		y = aux_y + vy*dt

		!Mudamos a densidade conforme a altitude
		r = ro*((1 - (y*(beta/To)))**alpha)

		aux = vx

		!Verificamos se a velocidade em x ou y está acima da velocidade do som e usamos o coef. de arrasto adequado

		IF ((vx < v_som) .AND. (vy < v_som)) THEN
	
			vx = vx - (((r*A*CD_menor)/(2*m))*sqrt((vx**2)+(vy**2))*vx*dt)
		
			vy = vy - ((g+(((r*A*CD_menor)/(2*m))*sqrt((aux**2)+(vy**2))*vy))*dt)
		ELSE IF ((vx > v_som) .AND. (vy < v_som)) THEN

			vx = vx - (((r*A*CD_maior)/(2*m))*sqrt((vx**2)+(vy**2))*vx*dt)
		
			vy = vy - ((g+(((r*A*CD_menor)/(2*m))*sqrt((aux**2)+(vy**2))*vy))*dt)
		ELSE IF ((vy > v_som) .AND. (vx < v_som)) THEN

			vx = vx - (((r*A*CD_menor)/(2*m))*sqrt((vx**2)+(vy**2))*vx*dt)
		
			vy = vy - ((g+(((r*A*CD_maior)/(2*m))*sqrt((aux**2)+(vy**2))*vy))*dt)

		ELSE

			vx = vx - (((r*A*CD_maior)/(2*m))*sqrt((vx**2)+(vy**2))*vx*dt)
		
			vy = vy - ((g+(((r*A*CD_maior)/(2*m))*sqrt((aux**2)+(vy**2))*vy))*dt)

		endif

		v = sqrt((vx*vx)+(vy*vy))
		t = t + dt

		write(37,*)t,v



	enddo

	write(36,*)x

	
	end program projetil_veloz
