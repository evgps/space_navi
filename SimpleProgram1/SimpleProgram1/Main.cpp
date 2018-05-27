// GSO1.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
/*Программа ориентирования КА на ГСО*/
#include "XKY_HBO_4_2.h"



double max_angle(KU_TimeDATA t, double p[6], double p_noisy[6], int tip)
{
	double DEL_UT1 = 0.; //- поправка в значение UTC
	double lon, alt;
	double lon_noisy, alt_noisy;
	double lon_err, alt_err;
	// пересчет в широту и долготу
	GEOCM(t, p, tip, DEL_UT1, &lon, &alt);
	GEOCM(t, p_noisy, tip, DEL_UT1, &lon_noisy, &alt_noisy);

	lon_err = abs((sin(lon - lon_noisy)));
	alt_err = abs(sin(alt - alt_noisy));
	printf("lon err: %lf\n", lon_err);
	printf("alt err: %lf\n", alt_err);
	return(fmax(lon_err, alt_err));
}

void make_consistent(double p[6])
{
	double PI = 3.141592653589793;
	p[1] = abs(p[1]);
	p[3] = abs(p[3]);
	for (int i = 2; i < 6; i++)
	{
		if (i == 3) continue;
		if (p[i] < 0)
			p[i] += 2 * PI;
		else if (p[i] > 2 * PI)
			p[i] -= 2 * PI;
	}



}


void add_randn(double p[6], double dP[6])
{
	int i,j;
	double e[6];
	int N_SAMPLES = 50;//ЦПТ в треде
	//Генерация ошибки ~ randn, используя центральную предельнуую теорему
	for (i = 0; i < 6; i++)
	{
		e[i] = 0;
		for (j = 0; j < N_SAMPLES; j++)
		{
			e[i] += 2.0*(rand()- RAND_MAX/2)/ RAND_MAX ;			//rand-0.5 ~ U[-1,1], mu=0, sigma = sqrt(4/12)
		}
		e[i] = (e[i] * sqrt(3. / N_SAMPLES)) * dP[i]; //получим randn ~ N[0, dR[i]]
		printf("error %lf\n", e[i]);
		p[i] += e[i];
	}
}

int main(int argc, _TCHAR* argv[])
{

/********************************************************************

Функция		:INT2000
Назначение	:Численное интегрирование уравнений движения КА методом
:Рунге-Кутты 4-го порядка

Возврат		: 0  - нормальное завершение работы функции
:-1  - начальный или конечный момент времени меньше нуля
:-2  - синус наклонения равен нулю
:-3  - радиус орбиты КА меньше экватор. радиуса Земли
:-4  - эксцентриситет орбиты КА больше 1

/*   typedef struct                                       */
/*       {                                                */
/*        int	d;                                        */
/*        double	s;                                    */
/*       }KU_TimeDATA;                                    */
/*                                          22.11.2005    */
/**********************************************************/

	int tip = 1;	// {Меридиальная СО (ГСО): 1,  экваториальная(ВЭО): 2}
	int strip = abs(tip);
	KU_TimeDATA ts; // начальный момент времени(сутки, секунды)
	KU_TimeDATA tt;	// конечный момент времени(сутки, секунды)

	//TimeDATA_u tt;	// конечный момент времени(сутки, секунды)
	int error = 0;

	int priz = 3; // переменная, задающая признаки учета возмущений
					// {прогнозирование кеплеровского движения КА : 1;
					// учет всех возмущений(потенциал, Луна, Солнце) : 3}
	int dot = 1;		// режим окончания процесса прогнозирования
					// { до заданного значения аргумента широты : 0;
					// до заданного значения : 1 }
	double p[6];
	double pn[6];	// начальный вектор состояний
	double pk[6];	// конечный вектор состояний
	double pn_noisy[6];	// начальный вектор состояний
	double pk_noisy[6];	// конечный вектор состояний

	double dR[6];
	double dP[6];

	double PI = 3.141592653589793;

	// Ссаная загрузка констант, без которой из-за FM=0 цикл в INTUM уходит в бесконечность, СПАСИБО ЧТО ПОСТАВИЛИ АССЕРТ

	coil();


	pn[0] = 42000.; //p - большая полуось  42 000 км
	pn[1] = 0.; //k = 0 (эксентриситет, а геостационрная орибат - это круговая орбита)  
	pn[2] = PI/4.; // q = pi/4 - это угол между направлениями из притягивающего центра на восходящий узел орбиты и на перицентр, мне кажется пох какой брать его
	pn[3] = 0.01; // i - наклонение не должно быть 0
	pn[4] = PI/6.; // omega - инерциальная долгота восходящего узла 
	pn[5] = PI/3.; //u

	// Проверка на геостационарность

	//Следует выбрать величину погрешности (сигму в N(0,сигма)
	// 1 matrix
	dR[0] = sqrt(1.675e+001);
	dR[1] = sqrt(3.620e-008);
	dR[2] = sqrt(1.509e-009);
	dR[3] = sqrt(7.088e-008);
	dR[4] = sqrt(7.246e-008);
	dR[5] = sqrt(1.418e-007);

	//// 2 matrix
	//dR[0] = sqrt(1.298e+001);
	//dR[1] = sqrt(3.309e-008);
	//dR[2] = sqrt(1.969e-009);
	//dR[3] = sqrt(1.060e-007);
	//dR[4] = sqrt(1.068e-007);
	//dR[5] = sqrt(1.210e-007);

	// 3 matrix
	//dR[0] = sqrt(2.375e+001);
	//dR[1] = sqrt(3.996e-008);
	//dR[2] = sqrt(1.237e-009);
	//dR[3] = sqrt(7.403e-008);
	//dR[4] = sqrt(7.302e-008);
	//dR[5] = sqrt(2.329e-007);


	//Конвертируем R в P
	RIPM(dR, dP);
	int i, t;
	double angle;
	int noisy_iter = 0;
	//Инициируем оба вектора
	for (int jj = 0; jj < 6; jj++)
		pn_noisy[jj] = pn[jj];

	// Добавим шумчику
	add_randn(pn_noisy, dR);
	ts.d = 0;
	ts.s = 0.;
	tt.d = 0;

	// Порог ошибки для угла между векторами (шумным и нешумным)
	double THRESHOLD = 0.1/180.*PI;
	printf("gogogo");
	double add_noise_interval = 28800.;
	double prop_s = 100.;
	for (int it = 0; it < 100000; it ++)
	{
		//интегрируем на prop_s секунды вперед
		tt.d = ts.d;
		tt.s = ts.s + prop_s;

		//каждые 8 часов добавляем шум
		if (it*(prop_s / add_noise_interval) > noisy_iter)
		{
			noisy_iter++;
			add_randn(pn_noisy, dR);
		}
		printf("Times: %d %lf %d %lf \n", ts.d, ts.s, tt.d, tt.s);
		error = INTUM(&ts, &tt, pn, pk, tip);
		if (error != 0) printf("First INTUM error: %d\n", error);
		
		error = INTUM(&ts, &tt, pn_noisy, pk_noisy, tip);
		if (error != 0)	printf("Second INTUM error: %d\n", error);

		//angle = calc_angle(pk, pk_noisy);

		angle = max_angle(tt, pk, pk_noisy, tip);
		
		//fprintf();
		if (abs(angle) > THRESHOLD)
		{
			printf("Angle is too hude: %lf \nEXIT!", angle);
			while (1);
			break;

		}
		else printf("Angle is smol and nice: %lf\n", angle);

		//копируем посчитанные pk в начальные
		for (int jj = 0; jj < 6; jj++)
		{
			pn[jj] = pk[jj];
			pn_noisy[jj] = pk_noisy[jj];
		}

		//переводим в дни
		ts.s += prop_s;
		ts.d += ts.s / 86400;
		ts.s = ts.s - 86400. *int(ts.s / 86400);

		
	}
	while (1);

}
