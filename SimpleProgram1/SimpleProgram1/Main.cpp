#include "stdafx.h"
/*Программа ориентирования КА на ГСО*/
#include "XKY_HBO_4_2.h"

// Очень грешно открываю файл глобально:
FILE * log_file;

// Поиск ошибок в пересчете на долготу и широту. Было лень нормально обрабатывать вычисления от 0 до 2pi поэтому asin(sin)
double max_angle(KU_TimeDATA t, double p[6], double p_noisy[6], int tip)
{
	double DEL_UT1 = 0.; //- поправка в значение UTC
	double lon, alt;
	double lon_noisy, alt_noisy;
	double lon_err, alt_err;
	// пересчет в широту и долготу
	GEOCM(t, p, tip, DEL_UT1, &lon, &alt);
	GEOCM(t, p_noisy, tip, DEL_UT1, &lon_noisy, &alt_noisy);

	lon_err = abs(asin(sin(lon - lon_noisy)));
	alt_err = abs(asin(sin(alt - alt_noisy)));
	printf("lon err: %lf\n", lon_err);
	printf("alt err: %lf\n", alt_err);
	fprintf(log_file, "%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", t.d, t.s, lon, alt, lon_noisy, alt_noisy, lon_err, alt_err, fmax(lon_err, alt_err));

	return(fmax(lon_err, alt_err));
}

// Проверка на говно
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

// Пошумим, блять!!!
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

// Страсть, 
// Cтраасть над ней имеет власть.
// Все громче звук,
// И понеслаааааась.
int main(int argc, _TCHAR* argv[])
{
	errno_t err;
	int tip = 1;	// {Меридиальная СО (ГСО): 1,  экваториальная(ВЭО): 2}
	int strip = abs(tip);
	KU_TimeDATA ts; // начальный момент времени(сутки, секунды)
	KU_TimeDATA tt;	// конечный момент времени(сутки, секунды)
	int error = 0;
	int priz = 3; // переменная, задающая признаки учета возмущений (кеплер - 1, все возмущ - 3)
	int dot = 1;		// режим окончания процесса прогнозирования (1 - до времени, 0 - до широты)
	double p[6];
	double pn[6];	// начальный вектор состояний
	double pk[6];	// конечный вектор состояний
	double pn_noisy[6];	// начальный вектор состояний
	double pk_noisy[6];	// конечный вектор состояний
	double dR[6]; //векторы ошибок
	double dP[6];
	double PI = 3.141592653589793;

	// Ссаная загрузка констант, без которой из-за FM=0 цикл в INTUM уходит в бесконечность, СПАСИБО ЧТО ПОСТАВИЛИ АССЕРТ
	coil();

////////////////////    ПЕРЕМЕННЫЕ ЭКСПЕРИМЕНТА  ///////////////////////////

	// Название файла
	err = fopen_s(&log_file, "experiment_0.csv", "a");

	// Начальные условия:
	pn[0] = 42000.; //p - большая полуось  42 000 км
	pn[1] = 0.; //k = 0 (эксентриситет, а геостационрная орибат - это круговая орбита)  
	pn[2] = PI/4.; // q = pi/4 - это угол между направлениями из притягивающего центра на восходящий узел орбиты и на перицентр, мне кажется пох какой брать его
	pn[3] = 0.01; // i - наклонение не должно быть 0
	pn[4] = PI/6.; // omega - инерциальная долгота восходящего узла 
	pn[5] = PI/3.; //u

	double THRESHOLD = 0.1 / 180.*PI;	// Порог ошибки для угла между векторами (шумным и нешумным)
	double add_noise_interval = 28800.;	// Интервал добавления шума (8ч)
	double prop_s = 100.;	// Шаг интегрирования ( интегрируем вперед на prop_s сек и сравниваем ошибки)

	//Следует выбрать величину погрешности (сигму в N(0,сигма))
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

	//// 3 matrix
	//dR[0] = sqrt(2.375e+001);
	//dR[1] = sqrt(3.996e-008);
	//dR[2] = sqrt(1.237e-009);
	//dR[3] = sqrt(7.403e-008);
	//dR[4] = sqrt(7.302e-008);
	//dR[5] = sqrt(2.329e-007);
///////////////////////////////////////////////////////////////////////////////////////////

	//Конвертируем R в P
	RIPM(dR, dP);

	fprintf(log_file, "p,k,q,i,omega,u\n%lf,%lf,%lf,%lf,%lf,%lf\n", pn[0], pn[1], pn[2], pn[3], pn[4], pn[5]);
	fprintf(log_file, "sigma_p,sigma_k,sigma_q,sigma_i,sigma_omega,u\n%lf,%lf,%lf,%lf,%lf,%lf\nt.d,t.s,lon,alt,lon_noisy,alt_noisy,lon_err,alt_err,fmax(lon_err, alt_err))\n", dP[0], dP[1], dP[2], dP[3], dP[4], dP[5]);

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

	for (int it = 0; it < 100500; it ++)
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
		printf("Time: %d %lf\n", ts.d, ts.s);

		// вычисления векторов p для 2х исходных векторов
		error = INTUM(&ts, &tt, pn, pk, tip);
		if (error != 0) printf("First INTUM error: %d\n", error);
		error = INTUM(&ts, &tt, pn_noisy, pk_noisy, tip);
		if (error != 0)	printf("Second INTUM error: %d\n", error);

		// нахождение ошибки
		angle = max_angle(tt, pk, pk_noisy, tip);
		
		// Выход при большом значении ошибки
		if (abs(angle) > THRESHOLD)
		{
			printf("Angle is too hude: %lf \nEXIT!", angle);
			break;
		}
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
	// закрыть все файлы
	int numclosed = _fcloseall();
	while (1);
}
