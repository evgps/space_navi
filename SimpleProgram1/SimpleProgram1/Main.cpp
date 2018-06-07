#include "stdafx.h"
#include "string.h"
/*Программа ориентирования КА на ГСО*/
#include "XKY_HBO_4_2.h"
#define _CRT_SECURE_NO_WARNINGS 1
typedef struct
{
	KU_TimeDATA	start;
	KU_TimeDATA	end;
	//stw[] - массив, состоящий из трех компонент :
	double stw[3];//s - управляющее ускорение по радиусу - вектору
	//t - управляющее ускорение по трансверсали
	//w - управляющее ускорение по бинормали

}Correction;

// Очень грешно открываю файл глобально:
FILE * log_file;
FILE * rand_file;
FILE * corr_file;
Correction * corrs;
// ОБОЖАЮ СИ ОХУЕННА СЧИТАТЬ СТРОКИ В ФАЙЛЕ ММММ
size_t count_lines(const char* filename) {
	FILE* fp;
	size_t cnt = 0u;
	fopen_s(&fp, filename, "r");
	while (!feof(fp)) {
		fscanf_s(fp, "%*[^\n]%*c");
		cnt++;
	}
	fclose(fp);
	return cnt;
}
// Секунд до
double secs_to(KU_TimeDATA from, KU_TimeDATA to)
{
	int days = to.d - from.d;
	double secs = days * 86400 + to.s - from.s;
	return secs;
}
// На вход KU_TimeDATA начала, длительность 
int read_corrections(const char * filename)
{
	int err;
	int n_str = count_lines(filename);
	printf("ANAL %d", n_str);
	// corr_file  is global variable 
	err = fopen_s(&corr_file, filename, "r");
	corrs = new Correction[n_str];
	
	//11715 52407.6139619252 11715 60822.6857554240 0.0000000000e+000 0.0000000000e+000 - 2.1766295620e-008
	for (int i = 0; i < n_str; i++)
	{
		if (fscanf_s(corr_file, "%d %lf %d %lf %le %le %le", &corrs[i].start.d, &corrs[i].start.s, \
			&corrs[i].end.d, &corrs[i].end.s,\
			&corrs[i].stw[0], &corrs[i].stw[1], &corrs[i].stw[2]) == NULL)
			printf("Deal with error.\n");
		printf("READ: %d %lf %d %lf %le %le %le\n", corrs[i].start.d, corrs[i].start.s, \
			corrs[i].end.d, corrs[i].end.s, \
			&corrs[i].stw[0], &corrs[i].stw[1], &corrs[i].stw[2]);
	}
	//*corrs = *tmp_corrs;
	err = fclose(corr_file);
	return(n_str);
}

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

// Пошумим, блять!!!
void add_randn(double p[6], double dP[6])
{
	int i,j;
	double e[6];
	int N_SAMPLES = 500;//ЦПТ в треде
	//Генерация ошибки ~ randn, используя центральную предельнуую теорему
	for (i = 0; i < 6; i++)
	{
		e[i] = 0;
		for (j = 0; j < N_SAMPLES; j++)
		{
			e[i] += 2.0*(rand()- RAND_MAX/2)/ RAND_MAX ;			//rand-0.5 ~ U[-1,1], mu=0, sigma = sqrt(4/12)
		}
		e[i] = (e[i] * sqrt(3. / N_SAMPLES)) * dP[i]; //получим randn ~ N[0, dR[i]]
		p[i] += e[i];
	}
}

int main(int argc, _TCHAR* argv[])
{
	errno_t err;
	int tip = 1;	// {Меридиальная СО (ГСО): 1,  экваториальная(ВЭО): 2}
	int strip = abs(tip);
	KU_TimeDATA ts; // начальный момент времени(сутки, секунды)
	KU_TimeDATA tt;	// конечный момент времени(сутки, секунды)
	KU_TimeDATA nextPOI;	// следующий важный момент (смена типа без корр/с корр)
	int INTU_type = 0; //1 - сейчас интервал с коррекцией. предположение о непересекаемости интервалов
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
	
	const char * filename = "PLAN_COR.TXT";
	// Название файла
	int cor_max_num;
	err = fopen_s(&log_file, "experiment_corr.csv", "a");
	cor_max_num = read_corrections(filename);

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
	double lat_stat, lon_stat;
	//Инициируем оба вектора
	for (int jj = 0; jj < 6; jj++)
		pk[jj] = pn[jj];
	for (int jj = 0; jj < 6; jj++)
		pn_init[jj] = pn[jj];

	ts.d = 11713;
	tt.d = ts.d;
	ts.s = 54880.;
	tt.s = ts.s - prop_s;

	GEOCM(tt, pn_init, tip, 0., &lat_stat, &lon_stat);
	printf("Стартовая позиция: %lf, %lf", lat_stat, lon_stat);
	nextPOI.d = corrs[0].start.d;
	nextPOI.s = corrs[0].start.s; // Считаем, что первое говно в файле коррекций позже старта!!!
	INTU_type = 0;			// Без коррекции 
	double secs_to_POI;
	double secs_integrate;
	int cor_num = 0; // Номер текущей коррекции
	TimeDATA_u tmp;
	for (int it = 0; it < 100500; it ++)
	{
		//интегрируем на secs_integrate секунды вперед
		tt.d = ts.d;
		
		secs_to_POI = secs_to(ts, nextPOI);
		secs_integrate = fmin(prop_s, secs_to_POI); // или квант или РТ (до нового события)
		printf("Time: %d %lf\n", ts.d, ts.s);

		tt.s = ts.s + secs_integrate;

		if (secs_integrate == 0)
		{
			if (INTU_type == 0) //переходим к коррекции ебала
			{
				printf("CHANGE INTU_type %d\n", INTU_type);

				nextPOI.d = corrs[cor_num].end.d;
				nextPOI.s = corrs[cor_num].end.s;
			}
			else
			{
				cor_num++;
				if (cor_num == cor_max_num - 1)
				{
					nextPOI.d = 300000; //бесконечна далека
				}
				else
				{
					printf("CHANGE NUM %d\n", cor_num);
					nextPOI.d = corrs[cor_num].start.d;
					nextPOI.s = corrs[cor_num].start.s;
				}

			}
			INTU_type++;
			INTU_type = INTU_type % 2;
			continue;
		}

		if (INTU_type == 0)
		{
			// вычисления векторов p для 2х исходных векторов
			error = INTUM(&ts, &tt, pn, pk, tip);
			if (error != 0) printf("First INTUM error: %d\n", error);
		}
		else
		{
			// вычисления векторов p для 2х исходных векторов
			tmp.d = tt.d;
			tmp.s = tt.s;
			tmp.u = 228; //ненужно
			printf("INTUS ANUS %d\n", cor_num);
			error = INTUS(&ts, &tmp, pk, 1, 1, tip, corrs[cor_num].stw);
			if (error != 0) printf("First INTUM error: %d\n", error);
		}

		// нахождение ошибки для нач точки
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
		}
		//переводим в дни
		ts.s += secs_integrate;
		ts.d += ts.s / 86400;
		ts.s = ts.s - 86400. *int(ts.s / 86400);
	}
	// закрыть все файлы
	int numclosed = _fcloseall();
	while (1);
}
