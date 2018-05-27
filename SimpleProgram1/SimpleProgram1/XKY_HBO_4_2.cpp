/***************************************************************************

                                XKY_HBO_4_2.c
                             -------------------
 
 ***************************************|************************************/
#include "stdafx.h"


//#include "XKY_HBO_5_0.h"
//#include "XKY_HBO_4_0.h"
//#include "XKY_HBO_4_1.h"
//#include "XKY_HBO_4_2.h"
//#include "XKY_HBO_4_6.h"
extern double MKC[];
extern double TKOC[];



/********************************************************************

	Функция		:INT2000
	Назначение	:Численное интегрирование уравнений движения КА методом 
                :Рунге-Кутты 4-го порядка
	Описание    :Исходные данные  находятся в файле INT2000.txt
	             Выходные данные  находятся в файле INT2000_out.txt

**-------------------------------------------------------------------
	Вход		:*t_n  - Указатель на стуктуру KU_TimeDATA (поле структуры
				:		 см. ниже), содержащую компоненты, задающие началь
				:		 ный момент времени (сутки,секунды).
				:		 интервала прогнозирования.
				:*t_k  - Указатель на стуктуру KU_TimeDATA,содержащую ком-
				:		 поненты, задающие конечный момент  времени,  если  
				:		 прогноз до заданного момента времени (dot = 1),			
				:        или указатель на стуктуру KU_TimeDATA, содержащую
				:		 компоненты, задающие количество витков и значение
				:		 аргумента  широты (dot = 0).
				: p[]  - вектор орбитальных параметров КА в начальной точке
	            :        интервала прогноза
				: dot  - переменная, определяющая режим окончания процесса
				:        прогнозирования
                :        dot = 1-признак прогноза до заданного значения
				:		 момента времени
                :        dot = 0-признак прогноза до заданного значения
				:		 аргумента широты
				:priz  - переменная, задающая признаки учета возмущений
				:   priz = 1 прогнозирование кеплеровского движения КА
				:   priz = 3 учет всех возмущений (потенциал,Луна,Солнце)
	            :tip   - переменная, определяющия систему координат
				:     tip = 1, то система координат меридианальная (ГСО)
				:     tip = 2, то система координат экваториальная (ВЭО)
	  Вход/Выход:	нет
	Выход		:p[]   - вектор орбитальных параметров КА в конечной точке
	            :         интервала прогноза
				:*t_k  - Указатель на стуктуру KU_TimeDATA, содержащую ком-
				:		 поненты конечного момента времени в результате про 
				:		 гноза до заданного аргумента широты (dot = 0)			
	Возврат		: 0  - нормальное завершение работы функции
	            :-1  - начальный или конечный момент времени меньше нуля    
	            :-2  - синус наклонения равен нулю
				:-3  - радиус орбиты КА меньше экватор. радиуса Земли
				:-4  - эксцентриситет орбиты КА больше 1       
	Ошибки		:	нет
   typedef struct                                              
       {                                                      
         int	d;                                              
        double	s;                                          
       }KU_TimeDATA;                                         
   typedef struct                                              
       {                                                      
         int	d;                                              
        double	s;                                          
        double	u;                                          
       }TimeDATA_u;                                         
*-------------------------------------------------------------------
	Вызовы		:	SMJ2000				fun_11_05.h
	            :   ZVEWS
	            :   UGOLM
	            :   TUZM
	            :   RKSHM
				:   SECUSM
	Переменные
	-глобальные	:	нет
	-статические:	нет
	-внешние	:	нет
	Макросы		:	нет
	Зависимости	:	нет
*******************************************************************/
int INT2000(KU_TimeDATA *t_n,TimeDATA_u *t_k, double p[],
			int dot,int priz,int tip)
{
 int sbr,nsig;
 int i;
 double pred[18];
 double a,b,w0,pp[10],uscu[3],mjd,ut;
 double RE,DVAPI,FM,DS;
 KU_TimeDATA tk,t1;
 KU_TimeDATA *ts;
 KU_TimeDATA t_s;

 //Проверьте, что вы вызвали coil() для загрузки констант
  RE    = TKOC[ 9];
  DVAPI = TKOC[11];
  FM    = TKOC[ 3];
  DS    = TKOC[52];
  /*  Контроль  входных данных */ 
 if ((t_n->d<0)            ||(t_n->s<0.0))         return -1;
 if ((dot==1)&&((t_k->d<0) ||(t_k->s  <0.0)))      return -1;
 if (sin(p[3])==0.0)                               return -2;
 if (p[0]/(1.0+p[2]*cos(p[5])+p[1]*sin(p[5]))<=RE)
 {
	 p[0] = p[0];
//	 return -3;
 }
 if (sqrt(p[1]*p[1]+p[2]*p[2])>=1.0)               return -4;

 ut = 0.0;
  /*  БЛОК ПРЕДВАРИТЕЛЬНЫХ ВЫЧИСЛЕНИЙ */

    t_s.d = t_n->d;
    t_s.s = t_n->s;
    ts    = & t_s;

    tk.d  = t_k->d;
    tk.s  = t_k->s;



 /* Вычисление признака направления прогноза
    (1 - прогноз вперед, -1 - прогноз назад */
    nsig = 0;
  if(dot==1) a  = (tk.d-ts->d)*86400.0+tk.s-ts->s;
  else       a  =  tk.d*DVAPI+tk.s-p[5];
  if (a>=0.0) nsig =  1;
  if (a<0.0) nsig = -1;


  sbr = 1;

  uscu[0]    =     p[5];  // u - аргумент широты
  uscu[1]    = sin(p[5]); // синус аргумента широты
  uscu[2]    = cos(p[5]); // косинус аргумента широты



  pred[ 0]  = sin(p[3]); // синус наклонения
  pred[ 1]  = cos(p[3]); // косинус наклонения
  pred[ 2]  = sin(p[4]); // синус долготы ВУО
  pred[ 3]  = cos(p[4]); // косинус долготы ВУО
  pred[ 4]  = p[0];      // фокальный параметр
  pred[ 5]  = p[1];      // k = e*sin(w) Компоненты
  pred[ 6]  = p[2];      // q = e*cos(w) вектора Лапласа

  pred[ 7]  = sqrt(FM/p[0])/p[0];  // w00

 /* Задание шага интегрирования в зависимости от типа орбиты 
                9 град - ВЭО   и  18 град - ГСО*/
   pred[ 8]  = DVAPI/((tip==1)?20.0:40.0)*nsig;

  /* Перевод системного времени в модиф.юлиан.дату и
      вычисление среднего звездного времени*/
  t1.d  = ts->d;
  t1.s  = 0.0;
  SMJ2000(&t1, &mjd);
  ZVEWS(mjd,0.0, &pred[9]);
//  ZVEWS(mjd,-0.51, &pred[9]);

  pred[10] = DS;        //0,017202791805307 Изменение Гринв.зв.вр. за сутки

  pred[15] = sqrt(a = p[1]*p[1]+p[2]*p[2]); // Вычисление эксцентриситета
  if (!pred[15]) w0=0.0;
  else UGOLM(p[1]/pred[15],p[2],&w0);   // Вычисление аргумента перигея
  pred[11] = DVAPI-w0;                  // Вычисление угла от перигея до ВУО КА
  if (pred[11]+uscu[0] >= DVAPI) pred[11]-=DVAPI;

  a        = p[0]/(1.0-a);             // Вычисление большой полуоси
  pred[16] = a*sqrt(a/FM);
  pred[12] = DVAPI*pred[16];           // Вычисление оскулирующего периода
  pred[13] = pred[14]=0.0;
  pred[17] = sqrt((1.0-pred[15])/(1.0+pred[15]));

  TUZM(0.0    ,&pred[13],pred);
  TUZM(uscu[0],&b       ,pred);
  pred[14] = ts->s-b;


 
  for (i=0;i<10;i++) pp[i]=0.0;

  /* Организация интегрирования по шагам методом Рунге-Кутты  */
  int iter = 0;
  while(1)    {
    if(dot==0) 
	{
		a=tk.d*DVAPI+tk.s-uscu[0];
	    if (fabs(a)<=fabs(pred[8]))  
		{ pred[8]=a;  break; }
	} 
    else 
	{ 
		a=(tk.d - ts->d)*86400.0 + tk.s;
		if((ts->s - a)*nsig >= 0.0)
		{ SECUSM(ts, a, &tk.d,  pp, uscu, tip, priz, pred);
		            break;  }
	} 
    ut= ut+pred[8];
	RKSHM (ts, &tk.d, pp, uscu, sbr, dot, tip, priz, pred);
	iter++;

  } 
    ut= ut+pred[8];
    RKSHM (ts, &tk.d,  pp, uscu, sbr, dot, tip, priz, pred);


  for (i=0;i<5;i++) p[i] += pp[i];
                    p[5]  = uscu[0];
    t_k->d=t_s.d;
    t_k->s=t_s.s;
    t_k->u=ut;
  return 0;
}
/************************************* end INT2000 ****************************************/

/**********************************************************/
/*                       INTUM                            */
/*   Процедура прогнозирования  до заданного времени      */
/*                                                        */
/*   Вход:  KU_TimeDATA*tn, KU_TimeDATA*tk,               */
/*          double pn[6], int tip ;                       */
/*                                                        */
/*   Выход: double pk[6] ;                                */
/*                                                        */
/*   typedef struct                                       */
/*       {                                                */
/*        int	d;                                        */
/*        double	s;                                    */
/*       }KU_TimeDATA;                                    */
/*                                          22.11.2005    */
/**********************************************************/
int INTUM(KU_TimeDATA *tn, KU_TimeDATA *tk, double pn[], 
		   double pk[], int tip)
{
  int stip, priz;
  KU_TimeDATA ts;
  TimeDATA_u tt;
  int dot;
  int i;
  double p[6];

  ts.d = tn -> d;
  ts.s = tn -> s;
  tt.d = tk -> d;
  tt.s = tk -> s;
  for(i = 0; i < 6; i++)  p[i] = pn[i];
  stip = abs(tip);
  dot = 1;
  priz = 3;

  i = INT2000(&ts, &tt,  p, dot, priz, stip);

  if (i < 0) return i;
  for(i = 0; i < 6; i++)  pk[i] = p[i];
 return 0;
}

/********************************* end INTUM *********************************************/


/**********************************************************/
/*                      INTUUM                            */
/*  Процедура прогнозирования до заданного аргумента      */
/*  широты                                                */
/*   Вход:  KU_TimeDATA*tn, int kvit,  double uk,         */
/*          double pn[6],   int tip;                      */
/*                                                        */
/*   Выход: KU_TimeDATA*tk, double pk[6] ;                */
/*                                                        */
/*   typedef struct                                       */
/*       {                                                */
/*        int	d;                                        */
/*        double	s;                                    */
/*       }KU_TimeDATA;                                    */
/*                                          22.11.2005    */
/**********************************************************/
int INTUUM(KU_TimeDATA *tn, KU_TimeDATA *tk, int kvit, double uk, 
		   double pn[],double pk[],int tip)
{
  int stip,priz;
  KU_TimeDATA ts;
  TimeDATA_u tt;
  int dot;
  int i;
  double p[6];

  ts.d = tn->d;
  ts.s = tn->s;
  tt.d = kvit;
  tt.s = uk;
 
  for(i = 0; i < 6; i++)  p[i] = pn[i];
  stip = abs(tip);
  dot = 0;
  priz = 3;
  i=INT2000(&ts, &tt,  p, dot, priz, stip);
  if (i < 0) return i;
  tk->d = tt.d;
  tk->s = tt.s;
   for(i = 0;i < 6; i++)  pk[i] = p[i];
 return 0;
}
/********************************* end INTUUM *********************************************/

/********************************************************************

	Функция		:   INTUS
	Назначение	:	Численное интегрирование уравнений движения КА методом 
                :   Рунге-Кутты 4-го порядка с учетом управляющих ускорений
	Описание    :	Исходные данные  находятся в файле INTUS.txt
	                Выходные данные  находятся в файле INTUS_out.txt

**-----------------------------------------------------------------------
	Вход		:	 *t_n  - Указатель на стуктуру KU_TimeDATA (поле структуры
				:			 см. ниже), содержащую компоненты, задающие началь
				:			 ный момент времени (сутки,секунды).
				:			 интервала прогнозирования.
				:   *t_k   - Указатель на стуктуру KU_TimeDATA, содержащую ком-
				:			 поненты, задающие конечный момент  времени,  если  
				:			 прогноз до заданного момента времени (dot = 1),			
				:            или указатель на стуктуру KU_TimeDATA, содержащую
				:			 компоненты, задающие количество витков и значение
				:			 аргумента  широты (dot = 0).
				:	 p[]   - вектор орбитальных параметров КА в начальной точке
	            :            интервала прогноза
				:    dot   - переменная, определяющая режим окончания процесса
				:            прогнозирования
                :            dot = 1-признак прогноза до заданного значения
				:			 момента времени
                :            dot = 0-признак прогноза до заданного значения
				:			 аргумента широты
				:    priz  - переменная, задающая признаки учета возмущений
				:            priz = 1 прогнозирование кеплеровского движения КА
				:            priz = 3 учет всех возмущений (потенциал, Луна, Солнце)
	            :     tip  - переменная, определяющия систему координат
				:            tip = 1, то система координат меридианальная (ГСО)
				:            tip = 2, то система координат экваториальная (ВЭО)
	            :     stw[]- массив, состоящий из трех компонент:
				:            s - управляющее ускорение по радиусу-вектору
				:            t - управляющее ускорение по трансверсали
				:            w - управляющее ускорение по бинормали
	  Вход/Выход:	нет
	Выход		:	p[]    - вектор орбитальных параметров КА в конечной точке
	            :            интервала прогноза
				:   *t_k   - Указатель на стуктуру KU_TimeDATA, содержащую ком-
				:			 поненты конечного момента времени в результате про 
				:			 гноза до заданного аргумента широты (dot = 0)			
	Возврат		:    0     - нормальное завершение работы функции
	            :   -1     - начальный или конечный момент времени меньше нуля    
	            :   -2     - синус наклонения равен нулю
				:   -3     - радиус орбиты КА меньше экваториального радиуса Земли
				:   -4     - эксцентриситет орбиты КА больше 1       
	Ошибки		:	нет
   typedef struct                                              
       {                                                      
         int	d;                                              
        double	s;                                          
       }KU_TimeDATA;                                         
   typedef struct                                              
       {                                                      
         int	d;                                              
        double	s;                                          
        double	u;                                          
       }TimeDATA_u;                                         
**-----------------------------------------------------------------------
	Вызовы		:	SMJ2000				fun_2000.h
	            :   ZVEWS
	            :   UGOLM
	            :   TUZM
	            :   RKSHUS
				:   SECUSS
	Переменные
	-глобальные	:	нет
	-статические:	нет
	-внешние	:	нет
	Макросы		:	нет
	Зависимости	:	нет
*************************************************************************/
int INTUS(KU_TimeDATA *t_n,TimeDATA_u *t_k, double p[],
			    int dot,int priz,int tip,double stw[3])
{
     int sbr,nsig,i;
     double pred[18],pp[10],uscu[6];
     double a,b,w0,mjd,ut;
     double RE,DVAPI,FM,DS;
     KU_TimeDATA tk,t1;
     KU_TimeDATA *ts;
     KU_TimeDATA t_s;

     RE    = TKOC[ 9];
     DVAPI = TKOC[11];
     FM    = TKOC[ 3];
     DS    = TKOC[52];


     /*  Контроль  входных данных  */
     if ((t_n->d<0)            ||(t_n->s<0.0))         return -1;
     if ((dot==1)&&((t_k->d<0) ||(t_k->s  <0.0)))      return -1;
     if (sin(p[3])==0.0)                               return -2;
   //  if (p[0]/(1.0+p[2]*cos(p[5])+p[1]*sin(p[5]))<=RE) return -3;
     if (sqrt(p[1]*p[1]+p[2]*p[2])>=1.0)               return -4;

     ut = 0.0;
 
     /*  БЛОК ПРЕДВАРИТЕЛЬНЫХ ВЫЧИСЛЕНИЙ */

	 t_s.d = t_n->d;
     t_s.s = t_n->s;
     ts    = & t_s;

     tk.d  = t_k->d;
     tk.s  = t_k->s;



     /* Вычисление признака направления прогноза
     (1 - прогноз вперед, -1 - прогноз назад */
     nsig = 0;
     if(dot==1)   a  = (tk.d-ts->d)*86400.0+tk.s-ts->s;
     else         a  =  tk.d*DVAPI+tk.s-p[5];
     if (a>0.0) nsig =  1;
     if (a<0.0) nsig = -1;


     sbr = 1;

     uscu[0]   =     p[5];  // u - аргумент широты
     uscu[1]   = sin(p[5]); // синус аргумента широты
     uscu[2]   = cos(p[5]); // косинус аргумента широты
     uscu[3]   = stw[0];
     uscu[4]   = stw[1];
     uscu[5]   = stw[2];


     pred[ 0]  = sin(p[3]); // синус наклонения
     pred[ 1]  = cos(p[3]); // косинус наклонения
     pred[ 2]  = sin(p[4]); // синус долготы ВУО
     pred[ 3]  = cos(p[4]); // косинус долготы ВУО
     pred[ 4]  = p[0];      // фокальный параметр
     pred[ 5]  = p[1];      // k = e*sin(w) Компоненты
     pred[ 6]  = p[2];      // q = e*cos(w) вектора Лапласа

     pred[ 7]  = sqrt(FM/p[0])/p[0];  // w00

     /* Задание шага интегрирования в зависимости от типа орбиты 
                  9 град - ВЭО   и  18 град - ГСО*/
     pred[ 8]  = DVAPI/((tip==1)?20.0:40.0)*nsig;

     /* Перевод системного времени в модиф.юлиан.дату и
        вычисление среднего звездного времени*/
     t1.d  = ts->d;
     t1.s  = 0.0;
     SMJ2000(&t1, &mjd);
     ZVEWS(mjd,0.0, &pred[9]);
//  ZVEWS(mjd,-0.51, &pred[9]);

     pred[10] = 1.7202791805307e-02;           // ds Изменение Гринв.зв.вр. за сутки

     pred[15] = sqrt(a = p[1]*p[1]+p[2]*p[2]); // Вычисление эксцентриситета
     if (!pred[15]) w0=0.0;
     else UGOLM(p[1]/pred[15],p[2],&w0);       // Вычисление аргумента перигея
     pred[11] = DVAPI-w0;                      // Вычисление угла от перигея до восх.узла орбиты КА
     if (pred[11]+uscu[0] >= DVAPI) pred[11]-=DVAPI;

     a        = p[0]/(1.0-a);                  // Вычисление большой полуоси
     pred[16] = a*sqrt(a/FM);
     pred[12] = DVAPI*pred[16];                // Вычисление оскулирующего периода
     pred[13] = pred[14]=0.0;
     pred[17] = sqrt((1.0-pred[15])/(1.0+pred[15]));

     TUZM(0.0    ,&pred[13],pred);
     TUZM(uscu[0],&b       ,pred);
     pred[14] = ts->s-b;


 
     for (i=0;i<10;i++) pp[i]=0.0;

     /* Организация интегрирования по шагам методом Рунге-Кутты 4-го пор.*/

      while(1)    {
        if(dot==0) 
		{
		  a=tk.d*DVAPI+tk.s-uscu[0];
	      if (fabs(a)<=fabs(pred[8]))  
		  { pred[8]=a;  break; }
		} 

        else 
		{ 
		  a=(tk.d - ts->d)*86400.0 + tk.s;
          if((ts->s - a)*nsig >= 0.0)
		  { SECUSS(ts, a, &tk.d,  pp, uscu, tip, priz, pred);
		                break;  }
		} 
        ut= ut+pred[8];
        RKSHUS(ts, &tk.d, pp, uscu, sbr, dot, tip, priz, pred);
	 } 

     ut= ut+pred[8];
     RKSHUS(ts, &tk.d,  pp, uscu, sbr, dot, tip, priz, pred);


     for (i=0;i<5;i++) p[i] += pp[i];
                       p[5]  = uscu[0];
     t_k->d=t_s.d;
     t_k->s=t_s.s;
     t_k->u=ut;
     return 0;
    }
/********************************** end INTUS *********************************************/

/**********************************************************/
/*                       INTUKOR                          */
/*   Процедура прогнозирования  до заданного времени      */
/*     с учетом управляющих ускорений                     */
/*   Вход:  KU_TimeDATA*tn, KU_TimeDATA*tk,               */
/*          double pn[6], int tip ,  double stw[3]        */
/*                                                        */
/*   Выход: double pk[6] ;                                */
/*                                                        */
/*   typedef struct                                       */
/*       {                                                */
/*        int	d;                                        */
/*        double	s;                                    */
/*       }KU_TimeDATA;                                    */
/*                                          06.09.2007    */
/**********************************************************/
int INTUKOR(KU_TimeDATA *tn, KU_TimeDATA *tk, double pn[], 
		   double pk[],int tip,double stw[3])
{
  int stip,priz;
  KU_TimeDATA ts;
  TimeDATA_u tt;
  int dot;
  int i;
  double p[6];

  ts.d=tn->d;
  ts.s=tn->s;
  tt.d=tk->d;
  tt.s=tk->s;
  for(i=0;i<6;i++)  p[i]=pn[i];
  stip  = abs(tip);
  dot  = 1;
  priz = 3;

  i=INTUS(&ts, &tt,  p, dot, priz, stip,stw);

  if (i<0) return i;
  for(i=0;i<6;i++)  pk[i]=p[i];
 return 0;
}
/********************************* end INTUKOR ********************************************/

/**********************************************************/
/*                      INTUUKOR                          */
/*  Процедура прогнозирования до заданного аргумента      */
/*  широты с учетом управляющих ускорений                 */
/*   Вход:  KU_TimeDATA*tn, int kvit, double uk,          */
/*          double pn[6],   int tip,  double stw[3]       */
/*                                                        */
/*   Выход: KU_TimeDATA*tk, double pk[6] ;                */
/*                                                        */
/*   typedef struct                                       */
/*       {                                                */
/*        int	d;                                        */
/*        double	s;                                    */
/*       }KU_TimeDATA;                                    */
/*                                          06.09.2007    */
/**********************************************************/
int INTUUKOR(KU_TimeDATA *tn, KU_TimeDATA *tk,
	          int kvit, double uk,double pn[],
			  double pk[],int tip,double stw[3])
{
  int stip,priz;
  KU_TimeDATA ts;
  TimeDATA_u tt;
  int dot;
  int i;
  double p[6];

  ts.d=tn->d;
  ts.s=tn->s;
  tt.d=kvit;
  tt.s=uk;
 
  for(i=0;i<6;i++)  p[i]=pn[i];
  stip  = abs(tip);
  dot  = 0;
  priz = 3;
  i=INTUS(&ts, &tt,  p, dot, priz, stip,stw);
  if (i<0) return i;
  tk->d = tt.d;
  tk->s = tt.s;
   for(i=0;i<6;i++)  pk[i]=p[i];
 return 0;
}
/********************************** end INTUUKOR ******************************************/

/*--------------------------------------------------------------------
|                               INTK                                 |
|  Процедура интегрирования на заданном интервале времени[tn, tk]    |
|   c учётом нескольких интервалов протяжённой коррекции             |
|                                                                    |
|   Вход:                                                            |
|          KU_TimeDATA* tn    - время начала интегрирования          |
|          KU_TimeDATA* tk    - время конца интегрирования           |
|          double       Pn[6] - вектор параметров орбиты в момент tn |
|          int          tip   - признак системы координат            |
|                               (1 - меридианальная с.к.             |
|                                2 - экваториальная с.к.)            |
|          KU_MKOR *ST_KOR    - структура массива коррекций          |
|                                                                    |
|   Выход:                                                           |
|          int          kik   - количество исполненных коррекций на  |
|                               интервале интегрирования             |
|          double       Pk[6] - вектор параметров орбиты в момент tk |
|                                                                    |
|   Обращается к проседурам : SCOVM  INTUKOR  INTK_1                 |
|                                                                    |
|                                                                    |
|   typedef struct                                                   |
|   {                                                                |
|       int	    d;       сутки                                       |
|       double  s;       секунды                                     |
|   }KU_TimeDATA;                                                    |
|                                                                    |
|   typedef struct                                                   |
|   {                                                                |
|     int         n_kor;          заданное число  коррекций          |
|     KU_TimeDATA Mtkor[360][2];  массив времён   коррекцй(tn и tk)  |
|     double      MW[360][3];     массив программ коррекцй(Ws Wt Wb) | 
|     int         TIPkor[360];    массив типа коррекции              |
|                               0 - нет коррекции, 1 - ODУ, 2 - DDУ  |
|   }KU_MKOR;         структура массива коррекций                    |
|                                            версия 02.06.2008       |
|-------------------------------------------------------------------*/
void INTK(KU_TimeDATA *tn, KU_TimeDATA *tk, double Pn[6], int tip, KU_MKOR *ST_KOR,
		      int *kik, double Pk[6])
{
	  int N;
      int Nkor, Nkor_max;

	  /* 02.06.2008 */
	  KU_TimeDATA MtKOR[360][2];
	  double MW[360][3];
      int    tip_kor[360];
      int    mi_kor[360];      /* массив индексов интервалов коррекции,
	                              попавших в [tn,tk]                   */

      KU_TimeDATA t_ndk_i,  t_kdk_i;
      KU_TimeDATA t_nkor, t_kkor;
      KU_TimeDATA ti_n,  ti_k;

      int i, j, l;
      int s1, s2, s3, s4, s5;

      double W[3];
      double W_0[3];
      double Pi[6];

	  /* 02.06.2008 */
      Nkor_max = 360;

	  /* распаковка структуры  ST_KOR */
	  N = ST_KOR->n_kor;
      for (j=0; j<N; j++)
	  {
		 MtKOR[j][0].d  = ST_KOR->Mtkor[j][0].d;
		 MtKOR[j][0].s  = ST_KOR->Mtkor[j][0].s;
		 MtKOR[j][1].d  = ST_KOR->Mtkor[j][1].d;
		 MtKOR[j][1].s  = ST_KOR->Mtkor[j][1].s;

		 MW[j][0]       = ST_KOR->MW[j][0];
		 MW[j][1]       = ST_KOR->MW[j][1];
		 MW[j][2]       = ST_KOR->MW[j][2];
         tip_kor[j]     = ST_KOR->TIPkor[j];
	  }
	  W_0[0] = 0.0;
	  W_0[1] = 0.0;
	  W_0[2] = 0.0;

      s1 = 0;
      s2 = 0;
      s3 = 0;
      s4 = 0;
      s5 = 0;

      *kik = 0;
      if (N == 0)              /* интервалов коррекции нет  */
	  {
          INTUKOR(tn, tk, Pn, Pk, tip, W_0);
          goto KONEC;
	  }

      i = 0;
	  Nkor = 0;
	  /* 02.06.2008 */
      for (j=0; j<Nkor_max; j++) mi_kor[j] = 0; 

     NC:
     /*  сравнение времён  */
	  t_ndk_i.d =  MtKOR[i][0].d;
	  t_ndk_i.s =  MtKOR[i][0].s;
	  t_kdk_i.d =  MtKOR[i][1].d;
	  t_kdk_i.s =  MtKOR[i][1].s;
 
      SCOVM (tn, tk,       &s1);      /* tn > tk      ->  s1 =  1 */
	                                  /* tn < tk      ->  s1 = -1 */
      SCOVM (tk, &t_kdk_i, &s2);      /* tk < tkdk_i  ->  s2 = -1 */
      SCOVM (tn, &t_ndk_i, &s3);      /* tn > tndk_i  ->  s3 =  1 */
      SCOVM (tn, &t_kdk_i, &s4);      /* tn < tkdk_i  ->  s4 = -1 */
      SCOVM (tk, &t_ndk_i, &s5);      /* tk > tndk_i  ->  s5 =  1 */

      if (((s1 == 1) && (s2 == -1) && (s3 == 1)) || ((s1 == -1) && (s4 == -1) && (s5 == 1)))
	  {
		 mi_kor[Nkor] = i;
         Nkor    = Nkor + 1;  /*  количество интервалов коррекции,
		                          попавших в [tn,tk]                */
	  }
      if (i < N -1)
	  {
          i = i + 1;
		  goto NC;
	  }
      *kik = Nkor;
	  
      if (Nkor == 0)    /* интервалов коррекции, попавших в [tn,tk] нет  */
	  {
          INTUKOR(tn, tk, Pn, Pk, tip, W_0);
          goto KONEC;
	  }
      if (Nkor == 1)
	  {
		 W[0] = MW[mi_kor[0]][0];
	     W[1] = MW[mi_kor[0]][1];
	     W[2] = MW[mi_kor[0]][2];
	     t_nkor.d =  MtKOR[mi_kor[0]][0].d;
	     t_nkor.s =  MtKOR[mi_kor[0]][0].s;
	     t_kkor.d =  MtKOR[mi_kor[0]][1].d;
	     t_kkor.s =  MtKOR[mi_kor[0]][1].s;

         INTK_1(tn, tk, Pn, tip, &t_nkor, &t_kkor, W, Pk);
         goto KONEC;
	  }

     /* начало цикла по интервалам коррекции, попавших в [tn,tk]  */
      l = 0;
     ML:
         if (s1 == -1)      /* tn < tk   */
		 {
	     t_nkor.d =  MtKOR[mi_kor[l]][0].d;
	     t_nkor.s =  MtKOR[mi_kor[l]][0].s;
	     t_kkor.d =  MtKOR[mi_kor[l]][1].d;
	     t_kkor.s =  MtKOR[mi_kor[l]][1].s;
		 }
		 else
		 {
	     t_nkor.d =  MtKOR[mi_kor[Nkor-l-1]][0].d;
	     t_nkor.s =  MtKOR[mi_kor[Nkor-l-1]][0].s;
	     t_kkor.d =  MtKOR[mi_kor[Nkor-l-1]][1].d;
	     t_kkor.s =  MtKOR[mi_kor[Nkor-l-1]][1].s;
		 }
      if (l == 0)           /* первый интервал */
	  {
         ti_n.d = tn->d;
         ti_n.s = tn->s;
         for (j=0; j<=5; j++) Pi[j] = Pn[j]; 
         if (s1 == -1)      /* tn < tk   */
		 {
	        ti_k.d =  MtKOR[mi_kor[l+1]][0].d;
	        ti_k.s =  MtKOR[mi_kor[l+1]][0].s;
		 }
		 else
		 {
	        ti_k.d =  MtKOR[mi_kor[Nkor-l-2]][1].d;
	        ti_k.s =  MtKOR[mi_kor[Nkor-l-2]][1].s;
		 }
         goto INT_1;
	  }
      if (l == Nkor-1)      /* последний интервал */
	  {
         ti_k.d = tk->d;
         ti_k.s = tk->s;
         if (s1 == -1)      /* tn < tk   */
		 {
	        ti_n.d =  MtKOR[mi_kor[l]][0].d;
	        ti_n.s =  MtKOR[mi_kor[l]][0].s;
		 }
		 else
		 {
	        ti_n.d =  MtKOR[mi_kor[Nkor-l-1]][1].d;
	        ti_n.s =  MtKOR[mi_kor[Nkor-l-1]][1].s;
		 }
         goto INT_1;
	  }
      else               /* внутренний интервал */
	  {
		 if (s1 == -1)      /* tn < tk   */
		 {
	        ti_n.d =  MtKOR[mi_kor[l]][0].d;
	        ti_n.s =  MtKOR[mi_kor[l]][0].s;
	        ti_k.d =  MtKOR[mi_kor[l+1]][0].d;
	        ti_k.s =  MtKOR[mi_kor[l+1]][0].s;
		 }
		 else
		 {
	        ti_n.d =  MtKOR[mi_kor[Nkor-l-1]][1].d;
	        ti_n.s =  MtKOR[mi_kor[Nkor-l-1]][1].s;
	        ti_k.d =  MtKOR[mi_kor[Nkor-l-2]][1].d;
	        ti_k.s =  MtKOR[mi_kor[Nkor-l-2]][1].s;
		 }
	  }
     INT_1:
         if (s1 == -1)      /* tn < tk   */
		 {
	  	    W[0] = MW[mi_kor[l]][0];
	        W[1] = MW[mi_kor[l]][1];
	        W[2] = MW[mi_kor[l]][2];
		 }
		 else
		 {
	  	    W[0] = MW[mi_kor[Nkor-l-1]][0];
	        W[1] = MW[mi_kor[Nkor-l-1]][1];
	        W[2] = MW[mi_kor[Nkor-l-1]][2];
		 }

         INTK_1(&ti_n, &ti_k, Pi, tip, &t_nkor, &t_kkor, W, Pk);
         for (j=0; j<=5; j++) Pi[j] = Pk[j];
		 
         if (l < Nkor-1)
		 {
		  l = l + 1;
          goto ML;
		 }
     KONEC: ;
}
/******************************* end INTK ************************************************/

/*--------------------------------------------------------------------
|                               INTK_1                               |
|   Процедура интегрирования на заданном интервале времени[tn, tk]   |
|        c учётом одного интервала протяжённой коррекции             | 
|                                                                    | 
|   Вход:                                                            | 
|      KU_TimeDATA* tn    - время начала интегрирования              | 
|      KU_TimeDATA* tk    - время конца интегрирования               | 
|      double       Pn[6] - вектор параметров орбиты в момент tn     | 
|      int          tip   - признак системы координат                | 
|                           (1 - меридианальная с.к.                 | 
|                            2 - экваториальная с.к.)                | 
|      KU_TimeDATA* tndk  - время начала действия  коррекции         | 
|      KU_TimeDATA* tkdk  - время конца действия  коррекции          |  
|      double       W[3]  - программа коррекции,действующая на борту |
|                                                                    |    
|   Выход:                                                           | 
|      double       Pk[6] - вектор параметров орбиты в момент tk     | 
|                                                                    | 
|   typedef struct                                                   | 
|    {                                                               |  
|        int	    d;        сутки                                  |          
|        double	s;        секунды                                    |      
|    }KU_TimeDATA;                                                   | 
|                                                                    |
|     Обращается к проседурам : SCOVM  INTUKOR                       |     
|                                                                    |
|                                               версия 27.07.2005    | 
|-------------------------------------------------------------------*/
void INTK_1(KU_TimeDATA *tn, KU_TimeDATA *tk, double Pn[6],int tip,
		       KU_TimeDATA *tndk, KU_TimeDATA *tkdk, double W[3],
		       double Pk[6])
{
    //  KU_TimeDATA t1,    t2;
      KU_TimeDATA t_kor_n,t_kor_k;
      double W_0[3];
	  double P1[6], P2[6];
      int s, s1, s2, s3, s4, s5;

	  W_0[0] = 0.0;
	  W_0[1] = 0.0;
	  W_0[2] = 0.0;

      SCOVM (tn, tk, &s);       /* tn < tk    ->  s  = -1 */
                                /* tn > tk    ->  s  =  1 */
      SCOVM (tk, tndk, &s1);    /* tk < tndk  ->  s1 = -1 */
                                /* tk = tndk  ->  s1 =  0 */
      SCOVM (tk, tkdk, &s2);    /* tk > tkdk  ->  s2 =  1 */
                                /* tk = tkdk  ->  s2 =  0 */
      SCOVM (tn, tndk, &s3);    /* tn < tndk  ->  s3 = -1 */
                                /* tn = tndk  ->  s3 =  0 */
      SCOVM (tn, tkdk, &s4);    /* tn > tkdk  ->  s4 =  1 */
                                /* tn = tkdk  ->  s4 =  0 */
      /* tn < tk */  
      if (s == -1)
	  {
         t_kor_n.d = tndk->d;
         t_kor_n.s = tndk->s;
         t_kor_k.d = tkdk->d;
         t_kor_k.s = tkdk->s;
	  }
	  else
	  {
         t_kor_n.d = tkdk->d;
         t_kor_n.s = tkdk->s;
         t_kor_k.d = tndk->d;
         t_kor_k.s = tndk->s;
	  }

      /* (tn < tk) и (tk <= tndk  или  tn >= tkdk) */
      if ((s == -1) && ((s1 == -1) ||(s1 == 0) || (s4 == 1) || (s4 == 0)))
	  {
         INTUKOR(tn, tk, Pn, Pk, tip, W_0);
         goto KONEC;
	  }
      /* (tn > tk) и (tn <= tndk  или tk >= tkdk) */
      if ((s == 1) && (((s3 == -1) ||(s3 == 0)) || ((s2 == 1) || (s2 == 0))))
	  {
         INTUKOR(tn, tk, Pn, Pk, tip, W_0);
         goto KONEC;
	  }

      SCOVM (tn, tk,   &s1);      /* tn > tk    ->  s1 =  1 */
	                              /* tn < tk    ->  s1 = -1 */
      SCOVM (tn, tndk, &s2);      /* tn < tndk  ->  s2 = -1 */
                                  /* tn = tndk  ->  s2 =  0 */
      SCOVM (tn, tkdk, &s3);      /* tn > tkdk  ->  s3 =  1 */
                                  /* tn = tkdk  ->  s3 =  0 */
      SCOVM (tk, tkdk, &s4);      /* tk < tkdk  ->  s4 = -1 */
                                  /* tk = tkdk  ->  s4 =  0 */
      SCOVM (tk, tndk, &s5);      /* tk > tndk  ->  s5 =  1 */
                                  /* tk = tndk  ->  s5 =  0 */
      /*  ((tn < tk ) и (tn <= tndk)) или  ((tn > tk ) и (tn >= tkdk))    */
      if  (((s == -1) && ((s2 ==  -1)||(s2 == 0))) || ((s == 1) && ((s3 == 1)||(s3 == 0))))
	  {
      /*  ((tn < tk)  и (tk <= tkdk)) или  ((tn > tk ) и (tk >= tndk))   */
      if  (((s == -1) && ((s4 ==  -1)||(s4 == 0))) || ((s == 1) && ((s5 == 1)||(s5 == 0))))    
	  {/* 1 n*/
            /*  ((tn < tk) и  (tn = tndk))  или ((tn > tk)  и (tn = tkdk))*/
		    if (((s == -1) && (s2 == 0)) || ((s == 1) && (s3 == 0)))
			{
              INTUKOR(tn,       tk, Pn, Pk, tip, W);
			}
			else
			{
              INTUKOR(tn,     &t_kor_n, Pn, P1, tip, W_0);
              INTUKOR(&t_kor_n, tk,     P1, Pk, tip, W);
			}
	  }/* 1 k*/ 
      else
	  { /* 2 n*/     
            /*  ((tn < tk) и  (tn = tndk))  или ((tn > tk)  и (tn = tkdk))*/
		    if (((s == -1) && (s2 == 0)) || ((s == 1) && (s3 == 0)))
			{
	            /*   ((tn < tk)  и (tk = tkdk)) или ((tn > tk) и (tk = tndk))*/	        
				if (((s == -1) && (s4 == 0)) || ((s == 1) && (s5 == 0)))
				{
				   INTUKOR(tn,       tk,       Pn, Pk, tip, W);
				}
				else
				{
                   INTUKOR(tn,       &t_kor_k, Pn, P2, tip, W);
                   INTUKOR(&t_kor_k,  tk,      P2, Pk, tip, W_0);
				}
			}
			else
			{
	            /*   ((tn < tk)  и (tk = tkdk)) или ((tn > tk) и (tk = tndk))*/	        
				if (((s == -1) && (s4 == 0)) || ((s == 1) && (s5 == 0)))
				{
					INTUKOR( tn,      &t_kor_n, Pn, P1, tip, W_0);			
					INTUKOR(&t_kor_n,  tk,      P1, Pk, tip, W);
				}
				else
				{
				   INTUKOR( tn,      &t_kor_n, Pn, P1, tip, W_0);			
                   INTUKOR(&t_kor_n, &t_kor_k, P1, P2, tip, W);
                   INTUKOR(&t_kor_k,  tk,      P2, Pk, tip, W_0);
				}
			}
	  }
	  }/* 2 k*/
	  else
	  {     
         /*  ((tn < tk)  и (tk <= tkdk)) или  ((tn > tk ) и (tk >= tndk))   */
         if  (((s == -1) && ((s4 ==  -1)||(s4 == 0))) || ((s == 1) && ((s5 == 1)||(s5 == 0))))
		 {/* 3 n*/
              INTUKOR(tn,       tk, Pn, Pk, tip, W);
		 }/* 3 k*/
	     else
		 {/* 4 n*/        
	        /*   ((tn < tk)  и (tk = tkdk)) или ((tn > tk) и (tk = tndk))*/	        
			if (((s == -1) && (s4 == 0)) || ((s == 1) && (s5 == 0)))
			{
				   INTUKOR(tn,       tk,       Pn, Pk, tip, W);
			}
			else
			{
                   INTUKOR(tn,       &t_kor_k, Pn, P2, tip, W);
                   INTUKOR(&t_kor_k,  tk,      P2, Pk, tip, W_0);
			}
		 }/* 4 k*/
	  }

     KONEC: ;
}
/****************************** end INTK_1 ***********************************************/

/**********************************************************************/
/*                               INTKU                                */
/*  Процедура прогнозирования до заданного аргумента широты [kv,uk]   */
/*  c учётом интервалов управляющих воздействий                       */
/*                                                                    */
/*   Вход:                                                            */
/*          KU_TimeDATA* t    - время начала прогнозирования          */
/*          double       P[6] - вектор параметров орбиты в момент t   */
/*          kv                - количество витков прогнозирования     */
/*          uk                - аргумент широты, до которого          */
/*                              необходимо спрогнозировать            */
/*          int          tip   - признак системы координат            */
/*                               (1 - меридианальная с.к.             */
/*                                2 - экваториальная с.к.)            */
/*          KU_MKOR *STK       - структура массива коррекций          */
/*                                                                    */
/*   Выход:                                                           */
/*          KU_TimeDATA* tk   - время окончания прогнозирования       */
/*          double      Pk[6] - вектор параметров орбиты в момент tk  */
/*                                                                    */
/*   Обращается к процедурам : SCOVM  INTUM INTUM INTUKOR  INTUUKOR   */
/*                                                                    */
/*                                                                    */
/*   typedef struct                                                   */
/*   {                                                                */
/*       int	    d;       сутки                                    */
/*       double  s;       секунды                                     */
/*   }KU_TimeDATA;                                                    */
/*                                                                    */
/*   typedef struct                                                   */
/*   {                                                                */
/*       int	 d;       сутки                                       */
/*       double  s;       секунды                                     */
/*       double  u;       аргумент широты                             */
/*   }TimeDATA_u;                                                     */
/*                                                                    */
/*   typedef struct                                                   */
/*   {                                                                */
/*     int n_kor;                 заданное число управляющих воздейст.*/
/*     KU_TimeDATA Mtkor[360][2]; массив времён управляющих           */
/*                                воздействий(tn и tk)                */
/*     double MW[360][3];         массив программ управляющих         */ 
/*                                воздействий(Ws Wt Wb)               */
/*     int   TIPkor[360];         массив типа коррекции               */
/*                                0 - нет коррекции, 1 - ODУ, 2 - DDУ */
/*   }KU_MKOR;         структура массива управляющих воздействий      */
/*                                            версия 02.06.2008       */
/**********************************************************************/
int INTKU(KU_TimeDATA *t ,double p[6],int kv,double uk,int tip,
		       KU_TimeDATA *tk,double pk[6],KU_MKOR *STK)
{
	  int i,j,n,sig,s1,s2,kor,nkor,pkor;
	  KU_TimeDATA tnk,tkk,tr;
	  TimeDATA_u tp;
	  double a,b,PI;
      double stw[3],pr[6],pz[6];

      PI = TKOC[1];
	  tip  = abs(tip);
      pkor = 0;
	  kor  = STK->n_kor;

      sig  = 0;
      a    = kv*2.0*PI + uk;
      b    = p[5];
	  if ((a-b)>0.0) sig =  1;
      if ((a-b)<0.0) sig = -1;

	  if ((kor==0)||(sig==0))
	  {
        tr.d = t->d; tr.s = t->s;
	    for (i=0; i<=5; i++)  pr[i]=p[i];
        tp.d = kv;
		goto M2;
	  }


      tnk.d = STK->Mtkor[0][0].d;
      tnk.s = STK->Mtkor[0][0].s;
      SCOVM (t,&tnk, &s1);            /* t > tnk      ->  s1 =  1 */
      tkk.d = STK->Mtkor[kor-1][1].d;
      tkk.s = STK->Mtkor[kor-1][1].s;
      SCOVM (t,&tkk, &s2);            /* t > tkk      ->  s2 =  1 */

	  if (((s1<=0)&&(sig<0))||((s2>=0)&&(sig>0)))
	  {
        tr.d = t->d; tr.s = t->s;
	    for (i=0; i<=5; i++)  pr[i]=p[i];
        tp.d = kv;
		goto M2;
	  }

      nkor = 0;
	  for (j=0; j<=kor-1; j++)
	  {
      tnk.d = STK->Mtkor[j][0].d;
      tnk.s = STK->Mtkor[j][0].s;
      SCOVM (t,&tnk, &s1);            /* t > tnk  ->  s1 =  1 */
      if (s1<0) {pkor=0; break;}
      tkk.d = STK->Mtkor[j][1].d;
      tkk.s = STK->Mtkor[j][1].s;
      SCOVM (t,&tkk, &s2);            /* t > tkk  ->  s2 =  1 */
      if (s2<0) {pkor=1; nkor++; break;}
	  nkor++;
	  }
      if (pkor==1&&sig>0) nkor--;
      tp.d = t->d; tp.s = t->s;
	  for (i=0; i<=5; i++)  pr[i]=p[i];

	  if (sig>0)
	  {
        for (j=nkor; j<=kor-1; j++)
		{
          tr.d = tp.d; tr.s = tp.s;
	      for (i=0; i<=5; i++)  pz[i]=pr[i];
		  if (pkor==0)
		  {
		   tp.d = STK->Mtkor[j][0].d;
           tp.s = STK->Mtkor[j][0].s;
           i=INT2000(&tr,&tp,pr,1,3,tip);
           if (i<0) return i;
           b = b + tp.u;
		   if(fabs(a-b)<1.0e-14) goto RET;
		   if((a-b)    < 0.0)
		   {
	        for (i=0; i<=5; i++)  pr[i]=pz[i];
            b = b - tp.u;  goto M1;
		   }
		  } 
          pkor=1;
		  tr.d = tp.d; tr.s = tp.s;
	      for (i=0; i<=5; i++)  pz[i]=pr[i];
		  stw[0] = STK->MW[j][0];
		  stw[1] = STK->MW[j][1];
		  stw[2] = STK->MW[j][2];
		  tp.d = STK->Mtkor[j][1].d;
          tp.s = STK->Mtkor[j][1].s;
          i=INTUS  (&tr,&tp,pr,1,3,tip,stw);
          if (i<0) return i;
          b = b + tp.u;
		  if(fabs(a-b)<1.0e-14) goto RET;
		  if((a-b)    < 0.0)
		  {
	       for (i=0; i<=5; i++)  pr[i]=pz[i];
           b = b - tp.u;  goto M1;
		  }
          pkor=0;
		}
	    tr.d = tp.d; tr.s = tp.s;
	  }

	  else
	  {
        for (j=nkor-1; j>=0; j--)
		{
          tr.d = tp.d; tr.s = tp.s;
	      for (i=0; i<=5; i++)  pz[i]=pr[i];
		  if (pkor==0)
		  {
		   tp.d = STK->Mtkor[j][1].d;
           tp.s = STK->Mtkor[j][1].s;
           i=INT2000(&tr,&tp,pr,1,3,tip);
           if (i<0) return i;
           b = b + tp.u;
		   if(fabs(a-b)<1.0e-14) goto RET;
		   if((a-b)    > 0.0)
		   {
	        for (i=0; i<=5; i++)  pr[i]=pz[i];
            b = b - tp.u;  goto M1;
		   }
		  }
          pkor=1;
		  tr.d = tp.d; tr.s = tp.s;
	      for (i=0; i<=5; i++)  pz[i]=pr[i];
		  stw[0] = STK->MW[j][0];
		  stw[1] = STK->MW[j][1];
		  stw[2] = STK->MW[j][2];
		  tp.d = STK->Mtkor[j][0].d;
          tp.s = STK->Mtkor[j][0].s;
          i=INTUS  (&tr,&tp,pr,1,3,tip,stw);
          if (i<0) return i;
          b = b + tp.u;
		  if(fabs(a-b)<1.0e-14) goto RET;
		  if((a-b)    > 0.0)
		  {
	       for (i=0; i<=5; i++)  pr[i]=pz[i];
           b = b - tp.u;  goto M1;
		  }
          pkor=0;
		}          
	    tr.d = tp.d; tr.s = tp.s;
	  }

  M1: n = kv-(int)floor(b/(2.0*PI));
      tp.d = n;

  M2: tp.s = uk;
      
	  if (pkor==0) 
	  {
        i=INT2000(&tr,&tp,pr,0,3,tip);
        if (i<0) return i;
	  }
      if (pkor==1)
	  {
        i=INTUS  (&tr,&tp,pr,0,3,tip,stw);
        if (i<0) return i;
	  }
 RET: tk->d = tp.d; tk->s = tp.s;
	  for (i=0; i<=5; i++)  pk[i]=pr[i];
      return 0;

}
/********************************** end INTKU ********************************************/

/*************************************************************/
/*                       APSIDM                              */
/*                                                           */
/*  Расчет вектора состояния в апсидальной точке,заданной на */
/*  любом витке до или после исходного вектора состояния КА. */
/*                                                           */
/*  При прогнозе в апогей АР=+1,а в перигей АР=-1,           */
/*  J- номер витка, на котором ищется апсидальная точка:     */
/* 	J>0 - вперед, J<0 - назад.                               */
/*  Начальный виток имеет J=0.                               */
/*  Витком вперед считается  траектория  между двумя         */
/*  последовательными узлами.                                */
/*                                                           */
/*  Вход :  tn      - время, на которое задан вектор         */
/*          pn[6]   - вектор состояния                       */
/*          tip     - тип орбиты                             */
/*          j       - виток                                  */
/*          AP      - 1 перигей, 1 апогей                    */
/*                                                           */
/*  Выход : tk      - время  в апогее/перигее                */
/*          pk[6]   - вектор в апогее/перигее                */ 
/*                                                           */
/*  Используются константы : PI, Eps_apsid, MaxIt            */ 
/*                                                           */
/*  Используются процедуры :  INTUUM()                       */
/*                                                           */
/*                                 16.11.2005                */
/*************************************************************/
int APSIDM(	/*Возвращаемое значение - признак сходимости (1 - сошлось)*/
			/*INPUT*/
			KU_TimeDATA* tn,/*сутки, на кот. задан вектор*/
			double pn[],/*вектор состояния*/
			int tip, /*( 1 стационар, 2 эллипс)*/
			int j,/*виток*/
			int ap,/* +1 - апогей, -1 - перигей*/
			/*OUTPUT*/
			KU_TimeDATA *tk,/*сутки в апогее/перигее*/
			double pk[]/*вектор в апогее/перигее*/
			)
{
	int i = 0;/*счетчик итераций*/
	int l,nvit,MaxIt;
	int Pr,k;
	double wpi, uk, e;
	double PI,DVAPI,Eps_apsid;
	double p[6];
	KU_TimeDATA t;

    PI       =TKOC[ 1];
    DVAPI    =TKOC[11];
	MaxIt    =(int)TKOC[49];
	Eps_apsid=TKOC[50];
	Pr = 1;  
	t.d = tn->d;
	t.s = tn->s;
	for (l = 0; l < 6; l++)  p[l] = pn[l];
    nvit=j;
	uk = 0.0;
	do
	{
	 i++; 
	 k=INTUUM( &t,tk,nvit,uk, p, pk,tip);
	 if(k<0) return k;
	 e = sqrt( pk[1] * pk[1] + pk[2] * pk[2] );
	 wpi = 0.;
	 if ( e>1.0e-8 )
		{
		wpi = acos(pk[2]/e);
		if (pk[1]<0.0) 	wpi = DVAPI - wpi;
		}
	 if (ap==-1)/*ищем  перигей*/ uk = wpi;
	 else 
	 {	if (ap==1)
			{/*апогей*/
			uk = wpi-PI;
			if (uk<0.0) uk += DVAPI;
			}
		else
		{
			printf ("APSIDM: illegal parametr ap\n");
			return(0);
		}
	 }
	 nvit = 0;
	 t.s = tk->s;
	 t.d = tk->d;
	 for (l=0; l<6; l++) p[l] = pk[l];
	}
	while (fabs(pk[5]-uk)>Eps_apsid&&i<MaxIt);
	if (i==MaxIt)  Pr = 0;
	return( Pr);
	}
/******************************* end APSIDM ***********************************************/

/*************************************************************/
/*                       APSIDK                              */
/*                                                           */
/*  Расчет вектора состояния в апсидальной точке,заданной на */
/*  любом витке до или после исходного вектора состояния КА. */
/*                                                           */
/*  При прогнозе в апогей АР=+1,а в перигей АР=-1,           */
/*  J- номер витка, на котором ищется апсидальная точка:     */
/* 	J>0 - вперед, J<0 - назад.                               */
/*  Начальный виток имеет J=0.                               */
/*  Витком вперед считается  траектория  между двумя         */
/*  последовательными узлами.                                */
/*                                                           */
/*  Вход :  tn      - время, на которое задан вектор         */
/*          pn[6]   - вектор состояния                       */
/*          tip     - тип орбиты                             */
/*          j       - виток                                  */
/*          AP      - -1 перигей, 1 апогей                    */
/*          STK     - структура  массива коррекций           */
/*                                                           */
/*  Выход : tk      - время  в апогее/перигее                */
/*          pk[6]   - вектор в апогее/перигее                */ 
/*                                                           */
/*  Используются константы : PI, Eps_apsid, MaxIt            */ 
/*                                                           */
/*  Используются процедуры :  INTKU()                        */
/*                                                           */
/*                                 02.06.2008                */
/*************************************************************/
int APSIDK(	/*Возвращаемое значение - признак сходимости (1 - сошлось)*/
			/*INPUT*/
			KU_TimeDATA* tn,/*сутки, на кот. задан вектор*/
			double pn[],/*вектор состояния*/
			int tip, /*( 1 стационар, 2 эллипс)*/
			int j,/*виток*/
			int ap,/* +1 - апогей, -1 - перигей*/
            KU_MKOR *STK,/* массив коррекций */
			/*OUTPUT*/
			KU_TimeDATA *tk,/*сутки в апогее/перигее*/
			double pk[]/*вектор в апогее/перигее*/
			)
{
	int i = 0;/*счетчик итераций*/
	int l,nvit,MaxIt;
	int Pr,k;
	double wpi, uk, e,DVAPI,PI,Eps_apsid;
	double p[6];
	KU_TimeDATA t;

       PI     = TKOC[ 1];
    DVAPI     = TKOC[11];
    MaxIt     = (int)TKOC[49];
    Eps_apsid = TKOC[50];
	Pr = 1;  
	t.d = tn->d;
	t.s = tn->s;
	for (l = 0; l < 6; l++)  p[l] = pn[l];
    nvit=j;
	uk = 0.0;
	do
	{
	 i++; 
     k=INTKU(&t,p,nvit,uk,tip,tk,pk,STK);
	 if(k<0) return k;
	 e = sqrt( pk[1] * pk[1] + pk[2] * pk[2] );
	 wpi = 0.;
	 if ( e > 1.0e-12 )
		{
		wpi = acos( pk[2] / e);
		if ( pk[1] < 0) 	wpi = DVAPI - wpi;
		}
	 if (ap == -1)/*ищем  перигей*/ uk = wpi;
	 else 
		if (ap == 1)
			{/*апогей*/
			uk = wpi - PI;
			if (uk < 0.) uk += DVAPI;
			}
		else
			{
		//	printf ("apsidm: illegal parametr ap\n");
			return(0);
			}
	 nvit = 0;
	 t.s = tk->s;
	 t.d = tk->d;
	 for (l=0; l<6; l++) p[l] = pk[l];
	}
	while (fabs( pk[5] - uk) > Eps_apsid && i < MaxIt);
	if (i == MaxIt)  Pr = 0;
	return( Pr);
}
/********************************* end APSIDK *********************************************/

/**********************************************************/
/*                       DRAKM                            */
/*   Расчет драконического периода                        */
/*                                                        */
/*   Вход : t, p[6], tip                                  */
/*   Выход: tdr                                           */
/*                                                        */
/*   Используемые процедуры :  INTUUM()                   */
/*                                       22.11.05         */
/**********************************************************/
int DRAKM(KU_TimeDATA t,double p[6],int tip,double *tdr)
{
     int i;
	  KU_TimeDATA t1,t2;
    double p1[6],p2[6];

	i=INTUUM(&t,&t1,0,0.0,p,p1,-tip);
    if (i<0) return i;
	i=INTUUM(&t,&t2,1,0.0,p,p2,-tip);
    if (i<0) return i;
    *tdr = (t2.d - t1.d)*86400.0 + (t2.s-t1.s);      
   if (i<0) return i;
   return 0;
}
/******************************* end  DRAKM ***********************************************/

/**********************************************************************/
/*                                 PROV                               */
/*          Программа, моделирующая частный алгоритм  PROV            */
/*            прогнозирования до заданного момента времени            */
/*                 с учетом управляющих воздействий                   */
/*  Вход:                                                             */
/*  int          tip   - признак системы координат                    */
/*  KU_TimeDATA  todn  - время привязки начального действ. вектора    */
/*  double      Podn[6]- начальный действующий вектор параметров КА   */
/*  KU_TimeDATA  t     - конечный момент времени                      */
/*  KU_MKOR  MDKOR -структура, содержащая массив действующих коррекций*/
/*  KU_MKOR  MPKOR -структура, содержащая массив планируемых коррекций*/
/*  Выход:                                                            */
/*  double       P[6]  - вектор параметров КА после прогноза  на      */
/*                       заданный момент времени                      */
/*   Обращается к процедурам :   INTK,   INTUM                        */
/*                                                                    */
/*                                            версия 01.04.2008       */
/**********************************************************************/
void PROV(KU_TimeDATA todn,double podn[6],int tip,KU_TimeDATA t,double p[6],
	        KU_MKOR  MDKOR,KU_MKOR  MPKOR)
{
	  int    i,kik,k_kord,k_korp;
	  double pkor[6];
	  KU_TimeDATA tkor;

       k_kord =  MDKOR.n_kor;
       k_korp =  MPKOR.n_kor;

       if((k_kord==0)&&(k_korp==0)) INTUM(&todn,&t,podn,p,tip);

       if (k_kord!=0)
	   {
	     tkor.d = MDKOR.Mtkor[k_kord-1][1].d;
		 tkor.s = MDKOR.Mtkor[k_kord-1][1].s;
	     SCOVM (&t,&tkor,&i);
         if (i<=0) INTK(&todn, &t,podn,tip,&MDKOR,&kik,p);
         else 
		 {
           INTK(&todn, &tkor,podn,tip,&MDKOR,&kik,pkor);
           if(k_korp!=0)INTK  (&tkor,&t,pkor,tip,&MPKOR,&kik,p);
           else          INTUM(&tkor,&t,pkor,p,tip);
		 }
       }
       else
	   {
	     if(k_korp!=0)
		 {
		   tkor.d = MPKOR.Mtkor[0][0].d;
		   tkor.s = MPKOR.Mtkor[0][0].s;
           SCOVM (&t,&tkor,&i);
           if (i==1) INTK (&todn,&t,podn,tip,&MDKOR,&kik,p);
           else      INTUM(&todn,&t,podn,p,tip);
		 }
	   }
   return;
}
/************************************ end PROV ********************************************/

/********************************************************************

	Функция		:  RKSHM
	Назначение	: Численное интегрирование уравнений движения КА методом 
                : Рунге-Кутты 4-го порядка на один шаг
	Описание    :	

    **-------------------------------------------------------------------
	Вход		: *t    - момент времени
				: dot   - переменная, определяющая режим окончания 
				:         процесса прогнозирования
				: nvit  - количество прогнозируемых витков КА 
	            :  tip  - переменная, определяющия систему координат
				: priz  - переменная, задающая признаки учета возмущений //тут -nan
  Вход/Выход    : нет
	Выход		:pp[]  -  вектор орбитальных параметров КА в конечной 
	            :         точке интервала прогноза
				: *t_k  - конечный момент времени интервала прогнозир.,
				:         если dot = 0
	Возврат		:	нет
	Ошибки		:	нет
    **-------------------------------------------------------------------
	Вызовы		:	PRAVM				fun_11_05.h
	            :   TUZM
	Переменные
	-глобальные	:	нет
	-статические:	нет
	-внешние	:	нет
	Макросы		:	нет
	Зависимости	:	нет
************************************************************************/
void RKSHM(KU_TimeDATA *t, int *nvit, double pp[], double uscu[],
	       int sbr, int dot, int tip, int priz,  double pred[])
{
  int i,j;
  double u0,aa[5],yy[10],rr[10],ab,ac,z[10],sec;
  double DVAPI;


  DVAPI = TKOC[11];
  u0    = uscu[0];
  aa[0] = aa[1]=aa[4]=pred[8]/2.0;
  aa[2] = aa[3]=pred[8];
  for (j=0;j<10;j++) yy[j]=rr[j]=pp[j];
  for (j=0;j<4;j++)
  {
	 PRAVM(t->d,rr,z,uscu, tip, priz, pred);
     ab=aa[j+1]/3.0;
     ac=aa[j];

     uscu[0]=u0+ac;
	 uscu[1]=sin(uscu[0]);
	 uscu[2]=cos(uscu[0]);

     for (i=0;i<10;i++)
	 {
	 pp[i] += ab*z[i];
	 rr[i]  = yy[i]+ac*z[i];
	 }
  }
  TUZM(uscu[0],&sec, pred);
  t->s = pp[5]+sec;

  if (sbr==1)
  { 
     if (uscu[0]<0.0)
	 {
	       uscu[ 0] += DVAPI;
	       pred[14] -= pred[12];
	       if (!dot) *nvit+=1;
	 }
     else  if (uscu[0]>=DVAPI)
	 {
		   uscu[ 0] -= DVAPI;
		   pred[14] += pred[12];
		   if (!dot) *nvit-=1;
	 }
     if (t->s<0.0)
	 {  
	       t->d-=1;
	       pred[14] += 86400.0;
	       t->s     += 86400.0;
	       pred[ 9] -= pred[10];
	 }
     else if (t->s>=86400.0)
	 {
	       t->d+=1;
	       pred[14] -= 86400.0;
	       t->s     -= 86400.0;
	       pred[ 9] += pred[10];
	 }
  }
}
/********************************************************************************************/

/********************************************************************************************/
 void SECUSM(KU_TimeDATA *t, double tk, int *data, 
	        double pp[],  double uscu[], int tip,
	        int priz, double pred[])
 {
   int i,j, sbr,dot;
   double epsh=1.e-10;
   double epsf=1.e-05;
   double x0,x,xh,su,cu,hh,f,f1,f2,f3,z,alfs,betas,pz[10];


     /* предварительные вычисления*/
   sbr = 0;
    dot = 1;
 	 x0 = uscu[0];
     su = uscu[1];
     cu = uscu[2];
    for(i=0; i<10; i++) pz[i]=pp[i];
     hh = pred[8];
     f3 = t->s-tk;
     xh = x0-hh;
      x = xh;

   for(j=1; j<=4; j++)
   {
     /* восстановление начальных данных для RKSHM, 
		вычисление шага интегрирования и обращение к RKSHM */
     uscu[1]=su; uscu[2]=cu;
     for(i=0; i<10; i++) pp[i]=pz[i];
	 pred[8]=x-x0;
     RKSHM (t,data,pp,uscu,sbr,dot, tip,priz, pred);
     uscu[0]=x0;
      /* вычисление функции z и контроль (значение функции
	     или интервала перемены знака) по заданной погрешности */
     z=t->s-tk;
     if(fabs(z)<=epsf||fabs(hh)<=epsh) goto M222;
     switch(j)
	 {
      /* формирование данных для проведения хорды*/
     case 1: f1=z;
	         alfs=hh/(1.0-f3/f1);
	         x=xh+alfs;
	         break;

     case 2: f2=z;
	         if( FSIGN(f1)!=FSIGN(f2) )
			 {
		      xh=xh+hh;
		      hh=-hh;
		      alfs=alfs+hh;
		      f=f1; f1=f3; f3=f;
			 }
	         f=1.0-f2/f1;
      /* контроль на возникновение особых ситуаций*/
	         if( fabs(f2)>=fabs(f1)||
		         fabs(alfs)>fabs(hh*f)||fabs(f)==0.0 )
			 {
      /* особая ситуация, формирование данных для хода
		 по методу деления интервала пополам */
		     xh=xh+alfs;
		     hh=hh-alfs;
		     f1=f2;
		     x=xh+hh/2.0;
			 }
		     else
			 {
	  /* формирование данных для проведения секущей*/
		      betas=alfs/f;
		      f1=f2;
		      x=xh+betas;
              j=3;
			 }
	         break;

	  /* сведение задачи с исходной после хода по 
		 методу деления интервала пополам*/
     case 3: hh=hh/2.0;
			 f=z;
	         if((f*f2)<=0.0) f3=f;
		     else
			 {
			  f1=f;
			  xh=x;
			 }
	         alfs=hh/(1.0-f3/f1);
	         x=xh+alfs;
	         j=1;
	         break;

	  /* сведение задачи с исходной после хода по 
		 методу секущей*/
     case 4: f2=z;
	         if( FSIGN(f1)==FSIGN(f2) )
			 {
		      f1=f2;
		      hh=hh-betas;
		      xh=xh+betas;
			 }
		     else
		     {
		      f3=f2;
		      hh=betas-alfs;
		      xh=xh   +alfs;
		     }
	         alfs=hh/(1.0-f3/f1);
	         x=xh+alfs;
	         j=1;
	         break;
     }
   }
     /* восстановление входных данных алгоритма*/
M222:
      uscu[1] = su;
      uscu[2] = cu;
		for( i=0; i<10; i++) pp[i]=pz[i];

 }

/********************************************************************************************/

/********************************************************************************************/
void TUZM(double u,double *sec,double pred[])
{
 double sig,fi,e;
 double PI,DVAPI,Eps_tuz;

  PI     =TKOC[ 1];
  DVAPI  = TKOC[11];
  Eps_tuz= TKOC[51];

  fi   = pred[11]+u;
  sig = floor(fi/DVAPI);
  if(fabs(fabs(fi)-DVAPI) <Eps_tuz)
  {      sig=FSIGN(fi);	  e=0.0;   } 
  else
  {  if(fabs(fabs(fi)-PI) <Eps_tuz)	e=fi;
     else  e=2.0*atan(pred[17]*tan(fi/2.0));
  }  
    if (e<0.0) e+=DVAPI;
    e-=pred[15]*sin(e);
  *sec=e*pred[16]+pred[12]*sig-pred[13]+pred[14];
}
/********************************************************************************************/

/********************************************************************

	Функция		:   RKSHUS
	Назначение	:	Численное интегрирование уравнений движения КА методом 
                :   Рунге-Кутты 4-го порядка на один шаг
	Описание    :	

    **-----------------------------------------------------------------------
	Вход		:	 *t   - указатель на стуктуру KU_TimeDATA, содержащую
	            :           момент времени в начале шага прогноза;
				:    dot  - переменная, определяющая режим окончания процесса
				:           прогнозирования
				:    nvit - количество прогнозируемых витков КА (используется
				:           в режиме прогнозирования до аргумента широты)
	   		    :	 pp[] - вектор интегрируемых приращений в начальной точке
	            :           шага прогноза
	   		    : usucu[] - массив из шести компонент: u - аргумент широты,
	            :           su,cu - синус и косинус аргумента широты,
				:           s,t,w - управляющие ускорения;
	   		    :  pred[] - массив параметров из блока предварительных вычисл:
				:   [0-1] - синус и косинус наклонения орбиты КА;
		        :   [2-3] - синус и косинус долготы ВУО КА;
				:   [4-6] - p,k,q;
				:   [  7] - множитель, используется функцией PRAVM;
				:   [  8] - шаг интегрирования;
				:   [  9] - среднее Гринвичское звездное время на 0 часов 
				:           Московского времени текущей даты;
				:   [ 10] - DS- изменение Гринвичского звездного времени за 
				:           сутки;
				:   [ 11] - истинная аномалия ВУО; 
				:   [ 12] - оскулирующий период; 
				:   [ 13] - интервал времени от перигея до ВУО; 
				:   [ 14] - момент времени прохождения КА ВУО; 
				:   [ 15] - эксцентриситет орбиты КА; 
				:   [ 16] - 1/n0, n0-среднее движение; 
				:   [ 17] - множитель, используется функцией TUZM;
	            :     tip - переменная, определяющия систему координат;
				:    priz - переменная, задающая признаки учета возмущений;
	  Вход/Выход:	нет
	Выход		:	 pp[] - вектор интегрируемых приращений в конечной точке
	            :           шага прогноза
	   		    : uscu[0] - u - аргумент широты;
				:pred[ 9] - среднее Гринвичское звездное время на 0 часов 
				:           Московского времени текущей даты;
				:pred[14] - момент времени прохождения КА ВУО; 
				:   *t    - указатель на стуктуру KU_TimeDATA, содержащую
				:           момент времени в конце шага прогноза;
				:   *nvit - количество прогнозируемых витков КА (используется
				:           в режиме прогнозирования до аргумента широты).
	Возврат		:	нет
	Ошибки		:	нет
    **-----------------------------------------------------------------------
	Вызовы		:	PRAVUS	fun_2000.h
	            :   TUZM
	Переменные
	-глобальные	:	нет
	-статические:	нет
	-внешние	:	нет
	Макросы		:	нет
	Зависимости	:	нет
*************************************************************************/
void RKSHUS(KU_TimeDATA *t, int *nvit, double pp[], double uscu[],
	           int sbr, int dot, int tip, int priz,  double pred[])
{
     int i,j;
     double u0,aa[5],yy[10],rr[10],ab,ac,z[10],sec;
     double DVAPI;


     DVAPI = TKOC[11];


     /* предварительные вычисления */
     u0    = uscu[0];
     aa[0] = aa[1]=aa[4]=pred[8]/2.0;
     aa[2] = aa[3]=pred[8];
     /* дублирование в рабочие массивы yy и rr вектора рр */
     for (j=0;j<10;j++) yy[j]=rr[j]=pp[j];

     /* организация цикла 4-х этапного метода Рунге-Кутты 4-го порядка*/
     for (j=0;j<4;j++)
	 {
	  PRAVUS(t->d,rr,z,uscu, tip, priz, pred);
      ab=aa[j+1]/3.0;
      ac=aa[j];

      uscu[0]=u0+ac;
	  uscu[1]=sin(uscu[0]);
	  uscu[2]=cos(uscu[0]);

      for (i=0;i<10;i++)
	  {
	   pp[i] += ab*z[i];
	   rr[i]  = yy[i]+ac*z[i];
	  }
	 }
     /* конец цикла*/
     /* вычисление времени после выполненого шага интегрирования */
     TUZM(uscu[0],&sec, pred);
     t->s = pp[5]+sec;
     /* анализ переменной sbr - разрешена "нормализация" времени и
       аргумента широты*/
     if (sbr==1)
	 { 
      if (uscu[0]<0.0)
	  {
	       uscu[ 0] += DVAPI;
	       pred[14] -= pred[12];
	       if (!dot) *nvit+=1;
	  }
      else  if (uscu[0]>=DVAPI)
	  {
		   uscu[ 0] -= DVAPI;
		   pred[14] += pred[12];
		   if (!dot) *nvit-=1;
	  }
      if (t->s<0.0)
	  {  
	       t->d-=1;
	       pred[14] += 86400.0;
	       t->s     += 86400.0;
	       pred[ 9] -= pred[10];
	  }
      else if (t->s>=86400.0)
	  {
	       t->d+=1;
	       pred[14] -= 86400.0;
	       t->s     -= 86400.0;
	       pred[ 9] += pred[10];
	  }
	 }
}
/*****************************************************************************************/

/********************************************************************

	Функция		:   SECUSS
	Назначение	:	Организация вычислений в режиме прогнозирования 
                :   до заданного момента времени
	Описание    :	

    **-----------------------------------------------------------------------
	Вход		:	 *t   - указатель на стуктуру KU_TimeDATA, задающую момент
	            :           времени в конце интервала перемены знака функции;
                :    tk   - конечный момент времени;
                :   *data - первая компонента заданного конечного момента
				:           времени;
	   		    :	 pp[] - вектор интегрируемых приращений 	            :           шага прогноза;
	   		    : usucu[] - массив из шести компонент: u - аргумент широты,
	            :           su,cu - синус и косинус аргумента широты,
				:           s,t,w - управляющие ускорения;
	            :     tip - переменная, определяющия систему координат;
				:    priz - переменная, задающая признаки учета возмущений;
	   		    :  pred[] - массив параметров из блока предварительных вычисл:
				:   [0-1] - синус и косинус наклонения орбиты КА;
		        :   [2-3] - синус и косинус долготы ВУО КА;
				:   [4-6] - p,k,q;
				:   [  7] - множитель, используется функцией PRAVM;
				:   [  8] - шаг интегрирования;
				:   [  9] - среднее Гринвичское звездное время на 0 часов 
				:           Московского времени текущей даты;
				:   [ 10] - DS- изменение Гринвичского звездного времени за 
				:           сутки;
				:   [ 11] - истинная аномалия ВУО; 
				:   [ 12] - оскулирующий период; 
				:   [ 13] - интервал времени от перигея до ВУО; 
				:   [ 14] - момент времени прохождения КА ВУО; 
				:   [ 15] - эксцентриситет орбиты КА; 
				:   [ 16] - 1/n0, n0-среднее движение; 
				:   [ 17] - множитель, используется функцией TUZM;
	  Вход/Выход:	нет
	Выход		:pred[ 8] = h - "усеченный" шаг по аргументу широты
	Возврат		:	нет
	Ошибки		:	нет
    **-----------------------------------------------------------------------
	Вызовы		:	RKSHUS				fun_2000.h
	Переменные
	-глобальные	:	нет
	-статические:	нет
	-внешние	:	нет
	Макросы		:	нет
	Зависимости	:	нет
*************************************************************************/
void SECUSS(KU_TimeDATA *t, double tk, int *data, 
	        double pp[],  double uscu[], int tip,
	        int priz, double pred[])
{
     int i,j, sbr,dot;
     double epsh=1.e-10;
     double epsf=1.e-05;
     double x0,x,xh,su,cu,hh,f,f1,f2,f3,z,alfs,betas,pz[10];

    /* предварительные вычисления*/
     sbr = 0;
     dot = 1;
 	 x0  = uscu[0];
     su  = uscu[1];
     cu  = uscu[2];
     for(i=0; i<10; i++) pz[i]=pp[i];
     hh  = pred[8];
     f3  = t->s-tk;
     xh  = x0-hh;
      x = xh;

     for(j=1; j<=4; j++)
     {
     /* восстановление начальных данных для RKSHM, 
		вычисление шага интегрирования и обращение к RKSHM */
	  uscu[1]=su; uscu[2]=cu;
      for(i=0; i<10; i++) pp[i]=pz[i];
	  pred[8]=x-x0;
	  RKSHUS(t,data,pp,uscu,sbr,dot, tip,priz, pred);
      uscu[0]=x0;
      /* вычисление функции z и контроль (значение функции
	     или интервала перемены знака) по заданной погрешности */
      z=t->s-tk;
      if(fabs(z)<=epsf||fabs(hh)<=epsh) goto M222;

      switch(j)
	  {
      /* формирование данных для проведения хорды*/
       case 1: f1=z;
	           alfs=hh/(1.0-f3/f1);
	           x=xh+alfs;
	           break;
 
       case 2: f2=z;
	           if( FSIGN(f1)!=FSIGN(f2) )
			   {
		        xh=xh+hh;
		        hh=-hh;
		        alfs=alfs+hh;
		        f=f1; f1=f3; f3=f;
			   }
	           f=1.0-f2/f1;
      /* контроль на возникновение особых ситуаций*/
			   if( fabs(f2)>=fabs(f1)||
		           fabs(alfs)>fabs(hh*f)||fabs(f)==0.0 )
			   {
      /* особая ситуация, формирование данных для хода
		 по методу деления интервала пополам */
				xh=xh+alfs;
		        hh=hh-alfs;
		        f1=f2;
		        x=xh+hh/2.0;
			   }
		       else
			   {
	  /* формирование данных для проведения секущей*/
		        betas=alfs/f;
		        f1=f2;
		        x=xh+betas;
                j=3;
			   }
	           break;

	  /* сведение задачи с исходной после хода по 
		 методу деления интервала пополам*/
       case 3: hh=hh/2.0;
			   f=z;
	           if((f*f2)<=0.0) f3=f;
		       else
			   {
			    f1=f;
			    xh=x;
			   }
	           alfs=hh/(1.0-f3/f1);
	           x=xh+alfs;
	           j=1;
	           break;

	  /* сведение задачи с исходной после хода по 
		 методу секущей*/
       case 4: f2=z;
	           if( FSIGN(f1)==FSIGN(f2) )
			   {
		        f1=f2;
		        hh=hh-betas;
		        xh=xh+betas;
			   }
		       else
			   {
		        f3=f2;
		        hh=betas-alfs;
		        xh=xh   +alfs;
			   }
	           alfs=hh/(1.0-f3/f1);
	           x=xh+alfs;
	           j=1;
	           break;
	  }
	 }

     /* восстановление входных данных алгоритма*/
     M222:
        uscu[1] = su;
        uscu[2] = cu;
	    for( i=0; i<10; i++) pp[i]=pz[i];
}
/*****************************************************************************************/

/*****************************************************************************************/
void PRAVM(int data,double rr[],double z[],double uscu[],
		   int tip, int   priz,double pred[])
{
 double si,ci,sw,cw,gs,gt,gw,sl,tl,wl,ss,ts,ws,ds,dt,dw,stw[10];
 double a,dp,dk,dq,sec,r1;
 double FM,fr;
 KU_TimeDATA t;

    FM = TKOC[ 3];
    fr = TKOC[58];
    dp = rr[0];
    dk = rr[1];
    dq = rr[2];

 rr[0] = pred[4] + dp;
 rr[1] = pred[5] + dk;
 rr[2] = pred[6] + dq;

     si = pred[0]+rr[6];
     ci = pred[1]+rr[7];
     sw = pred[2]+rr[8];
     cw = pred[3]+rr[9];

 stw[2] = cw*uscu[2]-sw*uscu[1]*ci;
 stw[0] = sw*uscu[2]+cw*uscu[1]*ci;
 stw[1] =    uscu[1]*si;

 stw[5] = si*sw;
 stw[3] =-si*cw;
 stw[4] = ci;

 stw[8] =-cw*uscu[1]-sw*uscu[2]*ci;
 stw[6] = cw*uscu[2]*ci-sw*uscu[1];
 stw[7] =    uscu[2]*si;

 if (tip == 2)
    {
	 fr = TKOC[59];
	 a    = stw[2];
     stw[2] = stw[1];
     stw[1] = stw[0];
     stw[0] = a;
	   a    = stw[5];
     stw[5] = stw[4];
     stw[4] = stw[3];
     stw[3] = a;
	   a    = stw[8];
     stw[8] = stw[7];
     stw[7] = stw[6];
     stw[6] = a;
    }
   r1     = 1.0+rr[1]*uscu[1]+rr[2]*uscu[2];
   stw[9] = rr[0]/r1;

   TUZM(uscu[0],&sec,pred);
   sec    = sec+rr[5];
   t.d = data;
   t.s = sec;

   gs=gt=gw=sl=tl=wl=ss=ts=ws=ds=dt=dw=0.0;

   if (priz==3) 
   {
    GEOM  (sec,&gs,&gt,&gw,pred[9],stw);
	LS2000(&t,0,&sl,&tl,&wl,stw);
    LS2000(&t,1,&ss,&ts,&ws,stw);
    SDSTW (&t,fr,  &ds,&dt,&dw,stw);
   }   

   a    = sqrt(rr[0]/FM);
   gs   = (gs+sl+ss+ds)*a;
   gt   = (gt+tl+ts+dt)*a;
   gw   = (gw+wl+ws+dw)*a;
   ws   = gw*uscu[1]*ci/(si*r1);
   ss   = 1.0/(r1*r1/(rr[0]*a)-ws);
   r1   = 1.0/r1;
   z[0] = 2.0*rr[0]*gt*r1*ss;
   z[1] = (((rr[1]+uscu[1])*r1+uscu[1])*gt-gs*uscu[2]-rr[2]*ws)*ss;
   z[2] = (((rr[2]+uscu[2])*r1+uscu[2])*gt+gs*uscu[1]+rr[1]*ws)*ss;
   z[3] = gw*uscu[2]*r1*ss;
   z[4] = gw*uscu[1]*r1*ss/si;
   z[6] = z[3]*ci;
   z[7] =-z[3]*si;
   z[8] = z[4]*cw;
   z[9] =-z[4]*sw;
   r1   = 1.0+pred[5]*uscu[1]+pred[6]*uscu[2];
   si   = pred[ 7]*r1*r1;
   dp  /= pred[4];
   ci   = (dk*uscu[1]+dq*uscu[2])/r1;
   sw   = ci*ci+2.0*ci;
   dp  *= dp*(1.875+dp*(-2.1875+dp*(2.4609375-2.70703125*dp)))-1.5;
   a    = sw + dp + sw*dp -              ws/si;
   z[5] =-a*(1.0-a*(1.0-a*(1.0-a*(1.0-a))))/si;
}
/******************************************************************************************/

/************************************************************************

	Функция		:   PRAVUS
	Назначение	:	Вычисление правых частей дифференциальных уравнений 
	Описание    :	

    **-----------------------------------------------------------------------
	Вход		:	data  - первая компонента системного времени
	   		    :	 rr[] - вектор интегрируемых приращений в начальной точке
	            :           шага прогноза
	   		    : usucu[] - массив из шести компонент: u - аргумент широты,
	            :           su,cu - синус и косинус аргумента широты,
				:           s,t,w - управляющие ускорения;
	   		    :  pred[] - массив параметров из блока предварительных вычисл:
				:   [0-1] - синус и косинус наклонения орбиты КА;
		        :   [2-3] - синус и косинус долготы ВУО КА;
				:   [4-6] - p,k,q;
				:   [  7] - множитель, используется функцией PRAVM;
				:   [  8] - шаг интегрирования;
				:   [  9] - среднее Гринвичское звездное время на 0 часов 
				:           Московского времени текущей даты;
				:   [ 10] - DS- изменение Гринвичского звездного времени за 
				:           сутки;
				:   [ 11] - истинная аномалия ВУО; 
				:   [ 12] - оскулирующий период; 
				:   [ 13] - интервал времени от перигея до ВУО; 
				:   [ 14] - момент времени прохождения КА ВУО; 
				:   [ 15] - эксцентриситет орбиты КА; 
				:   [ 16] - 1/n0, n0-среднее движение; 
				:   [ 17] - множитель, используется функцией TUZM;
	            :     tip - переменная, определяющия систему координат;
				:    priz - переменная, задающая признаки учета возмущений;
	  Вход/Выход:	нет
	Выход		:	  z[] - вектор орбитальных параметров КА в конечной точке
	            :           шага прогноза
	Возврат		:	нет
	Ошибки		:	нет
    **-----------------------------------------------------------------------
	Вызовы		:	TUZM				fun_11_05.h
	            :   GEOM
				:   LS2000
	Переменные
	-глобальные	:	нет
	-статические:	нет
	-внешние	:	нет
	Макросы		:	нет
	Зависимости	:	нет
*************************************************************************/
void PRAVUS(int data,double rr[],double z[],double uscu[],
		    int tip, int   priz,double pred[])
{
     double si,ci,sw,cw,gs,gt,gw,sl,tl,wl,ss,ts,ws,ds,dt,dw,stw[10];
     double a,dp,dk,dq,sec,r1,p,k,q;
     double FM,fr;
     KU_TimeDATA t;

     FM = TKOC[ 3];
     fr = TKOC[58];

     /* вычисление компонент p,k,q,si,ci,sw,cw, равных сумме соответствующих
	    компонент параметров ОКО и интегрируемых приращений */
     dp = rr[0];
     dk = rr[1];
     dq = rr[2];

     p  = pred[4] + dp;
     k  = pred[5] + dk;
     q  = pred[6] + dq;

     si = pred[0]+rr[6];
     ci = pred[1]+rr[7];
     sw = pred[2]+rr[8];
     cw = pred[3]+rr[9];

 /* вычисление элементов массива stw, имеющих
	 значение направляющих косинусов */
     stw[2] = cw*uscu[2]-sw*uscu[1]*ci;
     stw[0] = sw*uscu[2]+cw*uscu[1]*ci;
     stw[1] =    uscu[1]*si;

     stw[5] = si*sw;
     stw[3] =-si*cw;
     stw[4] = ci;

     stw[8] =-cw*uscu[1]-sw*uscu[2]*ci;
     stw[6] = cw*uscu[2]*ci-sw*uscu[1];
     stw[7] =    uscu[2]*si;

     if (tip == 2)
	 {
	  fr    = TKOC[59];
	   a    = stw[2];
     stw[2] = stw[1];
     stw[1] = stw[0];
     stw[0] = a;
	   a    = stw[5];
     stw[5] = stw[4];
     stw[4] = stw[3];
     stw[3] = a;
	   a    = stw[8];
     stw[8] = stw[7];
     stw[7] = stw[6];
     stw[6] = a;
	 }
     /* вычисление модуля радиуса-вектора КА*/
     r1     = 1.0+k*uscu[1]+q*uscu[2];
     stw[9] = p/r1;
     /* вычисление Московского времени, соответствующего значению
	    независимой переменной интегрирования */
     TUZM(uscu[0],&sec,pred);
     sec    = sec+rr[5];
     t.d = data;
     t.s = sec;

     gs=gt=gw=sl=tl=wl=ss=ts=ws=ds=dt=dw=0.0;

   if (priz==3) 
   {
    GEOM  (sec,&gs,&gt,&gw,pred[9],stw);
	LS2000(&t,0,&sl,&tl,&wl,stw);
    LS2000(&t,1,&ss,&ts,&ws,stw);
    SDSTW (&t, fr, &ds,&dt,&dw,stw);
   }   

	  /* нормирование ускорений*/
      a    = sqrt(p/FM);
      gs   = (gs+sl+ss+ds+uscu[3])*a;
      gt   = (gt+tl+ts+dt+uscu[4])*a;
      gw   = (gw+wl+ws+dw+uscu[5])*a;

     /* вычисление производной ss=dt/du*/
     ws   = gw*uscu[1]*ci/(si*r1);
     ss   = 1.0/(r1*r1/(p*a)-ws);

     /* вычисление вектора z, компонентами которого являются
	    производные интегрируемых параметров по аргументу широты*/
     r1   = 1.0/r1;
     z[0] = 2.0*p*gt*r1*ss;
     z[1] = (((k+uscu[1])*r1+uscu[1])*gt-gs*uscu[2]-q*ws)*ss;
     z[2] = (((q+uscu[2])*r1+uscu[2])*gt+gs*uscu[1]+k*ws)*ss;
     z[3] = gw*uscu[2]*r1*ss;
     z[4] = gw*uscu[1]*r1*ss/si;
     z[6] = z[3]*ci;
     z[7] =-z[3]*si;
     z[8] = z[4]*cw;
     z[9] =-z[4]*sw;

     /* вычисление компонеты z[5] - производной от времени*/
     r1   = 1.0+pred[5]*uscu[1]+pred[6]*uscu[2];
     si   = pred[ 7]*r1*r1;
     dp  /= pred[4];
     ci   = (dk*uscu[1]+dq*uscu[2])/r1;
     sw   = ci*ci+2.0*ci;
     dp  *= dp*(1.875+dp*(-2.1875+dp*(2.4609375-2.70703125*dp)))-1.5;
     a    = sw + dp + sw*dp -              ws/si;
     z[5] =-a*(1.0-a*(1.0-a*(1.0-a*(1.0-a))))/si;
}
/******************************************************************************************/

/******************************************************************************************/
void SDSTW(KU_TimeDATA *ts, double fr, double *s, double *t, double *w,
		    double stw[])
{
 int    i,psk,pr;
 double mjd,x[4],xs[4],xska[4],eka[4];
 double fk,gm,xi;
 double Rz;
  Rz    = TKOC[ 8];
 
  psk = 1;
  SMJ2000(ts, &mjd);
  SOLM (mjd,psk, x);
  for (i=0; i<=2; i++)
  {   eka[i] = stw[i]*stw[9]; 
       xs[i] =   x[i]*  x[3]; 
	 xska[i] =  xs[i]-eka[i];}
  
  xska[3] = sqrt(xska[0]*xska[0]+xska[1]*xska[1]+xska[2]*xska[2]);
  
  for (i=0; i<=2; i++) xska[i]=xska[i]/xska[3];

  
  pr = 1; 
  gm = (xs[0]*eka[0]+xs[1]*eka[1]+xs[2]*eka[2])/x[3];
  xi =  sqrt(stw[9]*stw[9]-Rz*Rz);
  if(gm<0.0) {
	  if(fabs(gm)>xi)  pr = 0;}

  fk = fr*pr/(xska[3]*xska[3]);

  *s = -(xska[0]*stw[0]+xska[1]*stw[1]+xska[2]*stw[2])*fk;
  *t = -(xska[0]*stw[6]+xska[1]*stw[7]+xska[2]*stw[8])*fk;
  *w = -(xska[0]*stw[3]+xska[1]*stw[4]+xska[2]*stw[5])*fk;
}
/******************************************************************************************/

/******************************************************************************************/
void LS2000(KU_TimeDATA *ts, int l, double *s, double *t, double *w,
		   double stw[])
{
 int    ind;
 double mu,mjd,x[4];
 double a,b,c,e,at;
 
  ind = 1;
  SMJ2000(ts, &mjd);
  if (l==0) {  mu=TKOC[5];       LUNM (mjd,ind, x); }
  if (l==1) {  mu=TKOC[4];       SOLM (mjd,ind, x); }

 a   = stw[9]/x[3];
 e   = x[0]*stw[0]+x[1]*stw[1]+x[2]*stw[2];
 b   = (a-2.0*e)*a;
 c   = b*(-1.5+b*(1.875  +b*(-2.1875+b*(2.4609375+
	   b*(-2.70703125    +b*(2.9326171875+
	   b*(-3.14208984375 +b*(3.338470458984-b*3.523941040039))))))));
 at  = mu/(x[3]*x[3]);
 b   = at*c;
 *s  = b*e-at*a*(1.0+c);
 *t  = b*(x[0]*stw[6]+x[1]*stw[7]+x[2]*stw[8]);
 *w  = b*(x[0]*stw[3]+x[1]*stw[4]+x[2]*stw[5]);
}
/******************************************************************************************/

/******************************************************************************************/
void GEOM(double sek, double *gs, double *gt, double *gw,
		  double sm0,double stw[])
{
 int    n,m,koe,nm;
 double sm,swt,cwt,sfi,cfi,c1,tgfi,sl,cl,ror,vrr,vlr,vfr;
 double pn2,dpn2,pn1,dpn1,pn,pmm,dpn,dpmm,skl,ckl,
	    frrn0,frn,e1,e2;
 double acr,adr,cckl,dskl,cskl,dckl,ssr,ccr,a;
 double OMEGAZ,RE,FM; 
 static double cn[]=
 { 0.0, 0.0, 0.0, 0.0,   -484164.95e-09,    0.05e-09, 2438.76e-09,
   957.16e-09, 2032.83e-09,  906.86e-09,  715.99e-09,  543.52e-09,
  -531.34e-09,  350.04e-09,  989.28e-09, -193.25e-09,   68.72e-09,
   -67.36e-09,  651.94e-09,	-451.45e-09, -293.19e-09,  191.57e-09,
  -143.45e-09,  -77.66e-09,   38.02e-09,   56.65e-09,  -90.89e-09,
  -268.71e-09,   14.97e-09,   90.39e-09,  280.07e-09,  321.10e-09,
   242.18e-09, -271.41e-09,    1.76e-09, -357.77e-09,   14.22e-09,
	48.83e-09,   23.55e-09,   69.74e-09,  -19.42e-09, -243.02e-09,
   -21.69e-09,  -66.39e-09,   68.57e-09, -113.68e-09,	28.62e-09,
   138.77e-09,   21.68e-09, -173.88e-09,   -7.05e-09,  -23.37e-09,
	69.32e-09, -120.84e-09,  177.44e-09,  -39.85e-09,   48.82e-09,
	71.69e-09,  -86.69e-09,   -2.55e-09,  -83.45e-09,  -42.27e-09,
   -35.77e-09,   11.50e-09,   45.02e-09,  127.02e-09,   97.25e-09,
   -50.35e-09,    8.81e-09,   17.64e-09,  -40.88e-09,  -39.90e-09,
	49.54e-09,    2.51e-09,   12.24e-09,   -5.63e-09,  -37.79e-09,
   -43.01e-09,   42.55e-09,   35.28e-09,  -52.12e-09,   11.57e-09,
    60.50e-09,  -67.42e-09,   40.00e-09,   12.46e-09,  -14.84e-09,
   -18.83e-09,   42.97e-09,   -2.30e-09,    9.86e-09,   -2.57e-09,
	43.74e-09,  -53.35e-09,   49.74e-09,  -26.59e-09,   -8.60e-09,
	63.13e-09,  -26.08e-09,    3.39e-09,   -6.60e-09,   19.15e-09, 
	36.27e-09,  -41.84e-09,  -31.76e-09,  -59.95e-09,  -22.44e-09,
   -14.15e-09,  -34.20e-09,   50.06e-09,    3.74e-09,   30.45e-09,
   -22.62e-09,   29.51e-09,  -35.21e-09,   24.45e-09,   33.55e-09,
	16.76e-09,    8.72e-09,   31.68e-09,  -52.74e-09,    0.52e-09,
	18.01e-09,  -19.21e-09,   46.30e-09,  -36.66e-09,   10.14e-09,	
    35.19e-09,   51.71e-09,  -36.57e-09,   10.05e-09,    6.64e-09,
	-1.02e-09,  -34.27e-09,  -28.53e-09,    8.88e-09,  -12.60e-09,
	-3.05e-09,   22.38e-09,  -16.99e-09,  -25.65e-09,   40.37e-09, 
	-4.84e-09,    2.78e-09,   -8.30e-09,  -29.02e-09,  -21.23e-09, 
   -17.25e-09,   18.75e-09,   20.33e-09,   12.97e-09,  -19.07e-09,
   -11.13e-09, -39.85e-09 };

 static double dn[]= 
 { 0.0, 0.0, 0.0, 0.0, 0.0, 0.01e-09, -1399.68e-09, 0.0,
   249.03e-09, -619.29e-09, 1400.74e-09,    0.0,      -471.95e-09,
   651.66e-09, -197.67e-09,  298.76e-09,    0.0,       -84.56e-09,
  -326.28e-09, -201.09e-09,   54.38e-09, -669.62e-09,    0.0,
    32.21e-09, -364.79e-09,    2.64e-09, -465.16e-09, -535.01e-09,
  -234.90e-09,    0.0,        86.66e-09,   96.87e-09, -244.93e-09,
  -121.81e-09,   15.35e-09,  150.66e-09,   32.62e-09,    0.0,
    58.52e-09,   52.70e-09,  -90.56e-09,   66.22e-09,   80.82e-09,
   309.82e-09,   73.74e-09,	 126.69e-09,    0.0,        18.98e-09,
   -27.27e-09,  -89.27e-09,   17.67e-09,  -63.18e-09,  217.16e-09,
   -97.63e-09,   -0.71e-09,	 105.63e-09,    0.0,      -124.29e-09,
   -34.03e-09, -152.56e-09,  -75.33e-09,  -50.28e-09,  -81.07e-09,
     0.11e-09,  -87.13e-09,  -41.91e-09,  -23.55e-09,    0.0,
   -35.76e-09, -104.05e-09, -151.73e-09,  -64.68e-09,   54.98e-09,
    34.53e-09,  -87.30e-09,   18.99e-09,   38.69e-09,  -27.92e-09,
   -50.30e-09,    0.0,       -35.65e-09,   23.17e-09,   28.97e-09,
     1.31e-09,   13.12e-09,   47.17e-09,   33.05e-09,   16.99e-09,
	23.80e-09,   33.02e-09,   -5.80e-09,  -11.40e-09,    0.0,
	41.50e-09,  -62.92e-09,   76.38e-09,  -10.97e-09,   62.72e-09,
	-2.74e-09,   -7.14e-09,   -4.46e-09,   46.96e-09,  -33.08e-09,
	 2.45e-09,   89.34e-09,   68.54e-09,    0.0,        28.61e-09,
	 0.84e-09,   15.27e-09,  -14.81e-09,  -17.91e-09,    7.74e-09,
	-5.08e-09,  -14.73e-09,   32.95e-09,    0.93e-09,  -36.70e-09,
   -32.70e-09,   45.16e-09,   -4.17e-09,    0.0,        10.65e-09,
   -31.15e-09,   10.58e-09,    3.67e-09,   11.32e-09,  -33.55e-09,
     1.04e-09,   26.09e-09,   41.10e-09,   19.61e-09,   20.31e-09,
	15.67e-09,   -4.51e-09,  -25.64e-09,  -11.36e-09,    0.0,
	36.22e-09,   35.47e-09,  -24.69e-09,   51.48e-09,    1.47e-09,
   -31.11e-09,  -11.26e-09,    5.70e-09,  -33.44e-09,   12.93e-09,
    -1.63e-09,    6.79e-09,    0.48e-09,  -39.08e-09,  -34.73e-09,
	7.85e-09 };	

  RE    = TKOC[ 9];
  OMEGAZ= TKOC[ 2];
  FM    = TKOC[ 3];

   sm  = OMEGAZ*sek+sm0;
   swt = sin(sm);
   cwt = cos(sm);
   sfi   = stw[2];
   cfi   = sqrt(fabs(1.0-sfi*sfi));
   c1  = 1.0/cfi;
   tgfi= sfi*c1;
   sl  = (stw[1]*cwt-stw[0]*swt)*c1;
   cl  = (stw[0]*cwt+stw[1]*swt)*c1;
   ror = RE/stw[9];
   vrr = 0.0;
   vlr = 0.0;
   vfr = 0.0;

/* группа операторов подготовки к началу рекуррентного
   процесса вычисления полиномов     */

    pn2 = 0.0;
   dpn2 = 0.0;
    pn1 = 1.0;
   dpn1 = 0.0;
     n  = 0;
	 m  = 0;
	 pn = sqrt(2.0);
	 pmm= pn;
	 dpn= 0.0;
	dpmm= 0.0;
	 skl= 0.0;
	 ckl= 1.0;
   frrn0= FM/(stw[9]*stw[9]);
   frn  = frrn0;

/* блок вычисления частных производных vrr,vfr,vlr от m-ой
   долготной гармоники */

A111:
   if(n!=m) 
   {
/* группа операторов вычисления нормированных функций Лежандра 
       и их производных  с индексами n и m */
	e1 = sqrt((4.0*n*n-1)/(n*n-m*m));
	e2 = sqrt((n*n-2.0*n-m*m+1)/(4.0*n*n-8.0*n+3));
	pn = e1*(sfi* pn1 - e2 * pn2      );
   dpn = e1*(sfi*dpn1 - e2 *dpn2 + pn1);
   frn = frn*ror;
   }
/*группа операторов вычисления частных производных потенциала
  vr/r  vfi/r  vl/r      */
   nm  = (n*(n+1))/2 + m + 1;
   acr = cn[nm];
   adr = dn[nm];

   cckl = acr*ckl;
   dskl = adr*skl;
   cskl = acr*skl;
   dckl = adr*ckl;

   vrr  = -frn*  (n+1)*(cckl+dskl)* pn     + vrr;
   vfr  =  frn*        (cckl+dskl)*dpn*cfi + vfr;
   vlr  = -frn* m *    (cskl-dckl)* pn     + vlr;

   koe = (n-16)*(m-16)*(n-m-16);
   if(koe!=0) 
   {
	   if(n!=m) 
	   { pn2 =  pn1;
	    dpn2 = dpn1;
		 pn1 =  pn;
		dpn1 = dpn;
	   }
	   n = n + 1; goto A111;
   }
   m = m + 1;
   if(m>16) 
   {
 /* Вычисление возмущающих ускорений */
   *gs = vrr;
   *gt = (vfr*cfi*stw[8] + vlr*stw[5])*c1*c1;
   *gw = (vfr*cfi*stw[5] - vlr*stw[8])*c1*c1;
   }
   else 
   {
      n = m;
	  a = sqrt((2.0*n+1.0)/(2.0*n));
	dpmm= a*(dpmm*cfi - pmm*tgfi);
	 pmm= a*pmm*cfi;
	 pn2= 0.0;
	dpn2= 0.0;
	  pn= pmm;
	 pn1= pmm;
	 dpn=dpmm;
	dpn1=dpmm;
   frrn0= frrn0*ror;
     frn= frrn0;
	 ssr= skl;
	 ccr= ckl;
	 skl= ssr*cl + ccr*sl;
	 ckl= ccr*cl - ssr*sl;
            goto A111;
   }
}
/********************************************************************************************/

