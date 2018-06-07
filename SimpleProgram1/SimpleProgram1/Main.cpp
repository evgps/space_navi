#include "stdafx.h"
/*Программа ориентирования КА на ГСО*/
#include "XKY_HBO_4_2.h"

// Очень грешно открываю файл глобально:
FILE * log_file;
FILE * rand_file;

// Поиск ошибок в пересчете на долготу и широту отн заданной точки
double max_angle_stat(KU_TimeDATA t, double p[6], double lat_stat, double lon_stat, int tip)
{
	double DEL_UT1 = 0.; //- поправка в значение UTC
	double lat, lon;
	double lat_err, lon_err;
	// пересчет в широту и долготу
	GEOCM(t, p, tip, DEL_UT1, &lat, &lon);

	lat_err = abs(asin(sin(lat - lat_stat)));
	lon_err = abs(asin(sin(lon - lon_stat)));
	printf("lat err: %lf\n", lat_err);
	printf("lon err: %lf\n", lon_err);
	fprintf(log_file, "%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", t.d, t.s, lat, lon, lat_stat, lon_stat, lat_err, lon_err, fmax(lat_err, lon_err));
	//fprintf(log_file, "%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", t.d, t.s, p[0], p[1], p[2], p[3], p[4], p[5]);

	return(fmax(lat_err, lon_err));
}

// Поиск ошибок в пересчете на долготу и широту. Было лень нормально обрабатывать вычисления от 0 до 2pi поэтому asin(sin)
double max_angle(KU_TimeDATA t, double p[6], double p_noisy[6], int tip)
{
	double DEL_UT1 = 0.; //- поправка в значение UTC
	double lat, lon;
	double lat_noisy, lon_noisy;
	double lat_err, lon_err;
	// пересчет в широту и долготу
	GEOCM(t, p, tip, DEL_UT1, &lat, &lon);
	GEOCM(t, p_noisy, tip, DEL_UT1, &lat_noisy, &lon_noisy);

	lat_err = abs(asin(sin(lat - lat_noisy)));
	lon_err = abs(asin(sin(lon - lon_noisy)));
	printf("lat err: %lf\n", lat_err);
	printf("lon err: %lf\n", lon_err);
	fprintf(log_file, "%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", t.d, t.s, lat, lon, lat_noisy, lon_noisy, lat_err, lon_err, fmax(lat_err, lon_err));
	//fprintf(log_file, "%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", t.d, t.s, p[0], p[1], p[2], p[3], p[4], p[5]);

	return(fmax(lat_err, lon_err));
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
	int N_SAMPLES = 500;//ЦПТ в треде

	//int err;
	//err=fopen_s(&rand_file, "rand_file_0.csv", "a");

	//printf("%d\n", err);
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
	//fprintf(rand_file, "%lf,%lf,%lf,%lf,%lf,%lf\n", e[0], e[1], e[2], e[3], e[4], e[5]);
	//int numclosed = _fcloseall();
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
	double pn_init[6];
	double dR[6]; //векторы ошибок
	double dP[6];
	double PI = 3.141592653589793;
	srand(304);
	// Ссаная загрузка констант, без которой из-за FM=0 цикл в INTUM уходит в бесконечность, СПАСИБО ЧТО ПОСТАВИЛИ АССЕРТ
	coil();

////////////////////    ПЕРЕМЕННЫЕ ЭКСПЕРИМЕНТА  ///////////////////////////

	// Название файла
	err = fopen_s(&log_file, "experiment_0.csv", "a");

	// Начальные условия:
	//pn[0] = 42000.; //p - большая полуось  42 000 км
	//pn[1] = 0.; //k = 0 (эксентриситет, а геостационрная орибат - это круговая орбита)  
	//pn[2] = PI/4.; // q = pi/4 - это угол между направлениями из притягивающего центра на восходящий узел орбиты и на перицентр, мне кажется пох какой брать его
	//pn[3] = 0.01; // i - наклонение не должно быть 0
	//pn[4] = PI/6.; // omega - инерциальная долгота восходящего узла 
	//pn[5] = PI/3.; //u



	pn[0] = 42165.5079267958;
	pn[1] = 0.0000144658;
	pn[2] = 0.0000336783;
	pn[3] = 1.5708271304;
	pn[4] = 1.5743659646;
	pn[5] = 0.7311609003;

	pn[0] = 42165.04066972509463;
	pn[1] = -0.00000022158530;
	pn[2] = 0.00003041071769;
	pn[3] = 1.57081023390739;
	pn[4] = 1.57124925859947;
	pn[5] = 0.47580873722654;

	double THRESHOLD = 0.1 * PI / 180.;	// Порог ошибки для угла между векторами (шумным и нешумным)
	double add_noise_interval = 28800.;	// Интервал добавления шума (8ч)
	double prop_s = 360;	// Шаг интегрирования ( интегрируем вперед на prop_s сек и сравниваем ошибки)

	//Следует выбрать величину погрешности (сигму в N(0,сигма))


	//// 2 matrix
	dR[0] = sqrt(1.298e+001);
	dR[1] = sqrt(3.309e-008);
	dR[2] = sqrt(1.969e-009);
	dR[3] = sqrt(1.060e-007);
	dR[4] = sqrt(1.068e-007);
	dR[5] = sqrt(1.210e-007);

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

	//fprintf(log_file, "p,k,q,i,omega,u\n%lf,%lf,%lf,%lf,%lf,%lf\n", pn[0], pn[1], pn[2], pn[3], pn[4], pn[5]);
	//fprintf(log_file, "sigma_p,sigma_k,sigma_q,sigma_i,sigma_omega,u\n%lf,%lf,%lf,%lf,%lf,%lf\nt.d,t.s,lat,lon,lat_noisy,lon_noisy,lat_err,lon_err,fmax(lat_err, lon_err))\n", dP[0], dP[1], dP[2], dP[3], dP[4], dP[5]);

	int i, t;
	double angle;
	int noisy_iter = 1;
	//Инициируем оба вектора
	for (int jj = 0; jj < 6; jj++)
		pn_noisy[jj] = pn[jj];

	for (int jj = 0; jj < 6; jj++)
		pn_init[jj] = pn[jj];

	double lat_stat, lon_stat;


	// Добавим шумчику
	add_randn(pn_noisy, dR);
	ts.d = 11710;
	tt.d = ts.d;
	ts.s = 46800.;
	tt.s = ts.s - prop_s;

	GEOCM(tt, pn_init, tip, 0., &lat_stat, &lon_stat);
	printf("Стартовая позиция: %lf, %lf", lat_stat, lon_stat);

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
		//angle = max_angle(ts, pk, pk_noisy, tip);
		//angle = max_angle(ts, pk, pn_init, tip);
		// Для нач точки
		angle = max_angle_stat(ts, pk, lat_stat, lon_stat, tip);
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
