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
FILE * p_file;
FILE * dR_file;
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
	
	for (int i = 0; i < n_str; i++)
	{
		if (fscanf_s(corr_file, "%d %lf %d %lf %le %le %le", &corrs[i].start.d, &corrs[i].start.s, \
			&corrs[i].end.d, &corrs[i].end.s, \
			&corrs[i].stw[0], &corrs[i].stw[1], &corrs[i].stw[2]) == EOF)
		{
			printf("Deal with error.\n");
			break;
		}
		printf("READ: %d %lf %d %lf %le %le %le\n", corrs[i].start.d, corrs[i].start.s, \
			corrs[i].end.d, corrs[i].end.s, \
			corrs[i].stw[0], corrs[i].stw[1], corrs[i].stw[2]);
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
	err = fopen_s(&log_file, "experiment_corr.csv", "w");
	cor_max_num = read_corrections(filename);

	// Начальные условия:
	fopen_s(&p_file, "vector.txt", "r");
	fscanf_s(p_file, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf", &pn[0], &pn[1], &pn[2], &pn[3], &pn[4], &pn[5]);
	fscanf_s(p_file, "%d\t%lf", &ts.d, &ts.s);
	fclose(p_file);
	tt.d = ts.d;
	tt.s = ts.s;
	double THRESHOLD = 0.10 * PI / 180.;	// Порог ошибки для угла между векторами (шумным и нешумным)
	double add_noise_interval = 28800.;	// Интервал добавления шума (8ч)
	double prop_s = 100;	// Шаг интегрирования ( интегрируем вперед на prop_s сек и сравниваем ошибки)

	double angle;
	double lat_stat, lon_stat;
	//Инициируем оба вектора
	for (int jj = 0; jj < 6; jj++)
		pk[jj] = pn[jj];
	for (int jj = 0; jj < 6; jj++)
		pn_init[jj] = pn[jj];



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
		tt.d += tt.s / 86400;
		tt.s = tt.s - 86400. *int(tt.s / 86400);

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
					nextPOI.d = 30000; //бесконечна далека
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
			error = INTUS(&ts, &tmp, pn, 1, 3, tip, corrs[cor_num].stw);
			//копируем посчитанные pk в начальные
			for (int jj = 0; jj < 6; jj++)
			{
				pk[jj] = pn[jj];
			}
			if (error != 0) printf("First INTUM error: %d\n", error);
		}

		// нахождение ошибки для нач точки
		angle = max_angle_stat(tt, pk, lat_stat, lon_stat, tip);
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
