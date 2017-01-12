#include <iostream>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <math.h>

using namespace std;

const int N=300;   // kuans_energyの分割の数(勝手に変更不可)
const int TRY=100000; // 試行回数
const double PI=3.141592;
const double DELTA=5.0; //ミラースリットの短辺の長さ(mm)
const double IPSI=60.0; //ミラースリットの長編の長さ(mm)
const double BETA=200; //遮蔽体の長さ(mm)
const double GAMMA=70.0; //ミラーの長さ(mm)
const double Z0=-(BETA+GAMMA); //z軸の0点の位置(mm)
const double NU=2000; //ミラーを出た後の長さ(mm)
const int PERBEAM=10; //実際の中性子源の１秒当たりの粒子数
const int MESTIME=100000; //測定時間(s)
const double MKG=1.67*1e-27; //中性子の質量(kg)
const double MMEV=939.6; //中性子の質量(Mev)
const double TOFM=1240; //TOF時の長さ(mm)
const double DT=0.0000001; //時間発展の刻み幅
const double XMIN=-0.5*IPSI; //スリット内にいる条件(壁に当たらない条件)
const double XMAX=0.5*IPSI;
const double YMIN=-0.5*DELTA;
const double YMAX=0.5*DELTA;
const double YZMAX=PI/30; //0.0128
const double XZMAX=PI/30; //0.0205
const double DALPHA=PI/180; //ミラーの傾き(z軸からｙ軸への角)
const double M1Y=YMAX;  //傾けた後のミラーのyz座標
const double M1Z=Z0+BETA;
const double M2Y=YMIN-DELTA*(1/cos(DALPHA)-1)+DELTA*tan(DALPHA)*sin(DALPHA);
const double M2Z=Z0+BETA+DELTA*sin(DALPHA);
const double M3Y=YMIN-DELTA*(1/cos(DALPHA)-1)+(DELTA*tan(DALPHA)+GAMMA)*sin(DALPHA);
const double M3Z=Z0+BETA+(DELTA*tan(DALPHA)+GAMMA)*cos(DALPHA);
const double M4Y=YMAX+GAMMA*sin(DALPHA);
const double M4Z=Z0+BETA+GAMMA*cos(DALPHA); 

void time_evolution(double* T, double* X,double*Y, double*Z, double VX,double VY,double VZ){
	*T=*T+DT;
	*X=*X+VX*DT;
	*Y=*Y+VY*DT;
	*Z=*Z+VZ*DT;
}

int main(){
    
	FILE* f=fopen("mirror_simulation4.txt","w");
	ifstream ifs("kuans_energy.txt"); // kuans_energyの読み込み
	
	srand((unsigned int)time(NULL));
	
	
	
	
	int k[N];         //必要ないが元ファイルkuans_energy.txtの左側を消すのが面倒なのでこれに入れておく
	int a;       //E軸方向の乱数(0~299)
	double b[N];      //kuans_energyの分布
	int tof;     //あるiでのTOF(0~3000mus)
	int d;       //counts軸の乱数(0~399)
	double azi;    //極角
	double pol;    //方位角
	double x, y, z;  //座標
	double vx,vy,vz; //速度成分
	int beamnumber=0; //粒子数
	int time=0; //時間
	double v; //中性子の速度(mm/s)
	double lambda;  //中性子の波長(オングストローム)
	double ene; //中性子のエネルギー(MeV)
	double t=0; //時間発展にかかった時間
	int ms=0; //ミラースリットに入った数
	int ref=0; //反射された中性子の数
	int s=0; //スリットを抜けた数
	double sine; //角度乱数生成時の重みづけに必要な変数
	double vh; //速度の水平成分(YZ成分)
	double dbeta; //最初の粒子の入射方向(z軸からｙ軸への角)
	
	for(int i=0; i<N; i++){
		ifs >> k[i] >> b[i];
	}
	
	for(int i=0; i<TRY; i++){
		
		int msok=0; //ミラースリットに入ったかどうか
		int rok=0; //反射したかどうか
		int sok=0; //ミラーを抜けたか
		
// energy
		
		for(int j=0; j<100000; j++){
			a=rand()%N;
			tof=(a+1)*(3000/N);
			d=rand()%400;		
			if(d<b[a]){break;}
			}
		
// degree
		
		for(int j=0; j<1000000; j++){
			sine=(double(rand())/RAND_MAX)*2-1;
			pol=asin(sine)+PI/2;
			//	pol=(double(rand())/RAND_MAX)*(0.0230)+(PI/2-0.0115);  //pi/2-0.0115~pi/2-0.0115  
			azi=(double(rand())/RAND_MAX)*(PI/16)-(PI/32);  //-pi/32~pi/32
			v=TOFM/(tof*1e-6);
			vx=v*sin(pol)*sin(azi);
			vy=v*cos(pol);
			vz=v*sin(pol)*cos(azi);
			if(-YZMAX<vy/vz && vy/vz<YZMAX && -XZMAX<vx/vz && vx/vz<XZMAX){break;}
			if(j==999999){cout << "number little" << endl;}
		}
		
//		cout << "pol" << pol  <<endl;
//		cout << "azi[" << i << "]=" << azi  <<endl;
		
// position
		x=(double(rand())/RAND_MAX)*IPSI-0.5*IPSI;
		y=(double(rand())/RAND_MAX)*DELTA-0.5*DELTA;
		z=Z0;
//		cout << "x[" << i << "]=" << x  <<endl;
//		cout << "y[" << i << "]=" << y  <<endl;

		beamnumber++;
		
		if(beamnumber%PERBEAM==0){time++;}
		if(time==MESTIME){break;}
			
		ene=0.5*MKG*v*v*1e-6*6.24*1e+18*1e-6;
		lambda=(2*PI*197.3)/(sqrt(2*MMEV*ene))*1e-5;
//		cout << "v[" << i << "]=" << v << endl;
//		cout << "ene[" << i << "]=" << ene << endl;
//		cout << "lambda[" << i << "]=" << lambda << endl;

// 遮蔽体の内部
		while(z<Z0+BETA && XMIN<x && x<XMAX && YMIN<y && y<YMAX){
				
			time_evolution(&t, &x, &y, &z ,vx ,vy ,vz);

			}	
			
		if(z>Z0+BETA){msok=1; ms++;}
//		cout << "x[" << i << "]=" << x << endl;
 //   	cout << "y[" << i << "]=" << y << endl;
 //     cout << "z[" << i << "]=" << z << endl;
//		cout << "vdt=" << v*DT << endl;
//		cout << "t[" << i << "]=" << t << endl;
//		cout << "msok[" << i << "]=" << msok << endl;
			
// ミラーの内部			
		if(msok==1){
			dbeta=atan(vy/vz);
			vh=sqrt(v*v-vx*vx);
			while(XMIN<x && x<XMAX && y<(((M1Y-M4Y)/(M1Z-M4Z))*(z-M1Z)+M1Y) && y<(((M3Y-M4Y)/(M3Z-M4Z))*(z-M3Z)+M3Y) && y>(((M2Y-M3Y)/(M2Z-M3Z))*(z-M2Z)+M2Y)){
					
				time_evolution(&t, &x, &y, &z ,vx ,vy ,vz);
				if(y<(((M2Y-M3Y)/(M2Z-M3Z))*(z-M2Z)+M2Y) && lambda>3 && 0<DALPHA-dbeta && DALPHA-dbeta<PI/180){
					rok=1; ref++; goto EXIT;  //反射
				}
			}
			EXIT: ;
				
			if(rok==1){
				vy=vh*sin(2*DALPHA-dbeta);
				vz=vh*cos(2*DALPHA-dbeta);
				while(XMIN<x && x<XMAX && y<(((M1Y-M4Y)/(M1Z-M4Z))*(z-M1Z)+M1Y) && y<(((M3Y-M4Y)/(M3Z-M4Z))*(z-M3Z)+M3Y)){
					time_evolution(&t, &x, &y, &z ,vx ,vy ,vz);
				}
			}
			if(y>(((M3Y-M4Y)/(M3Z-M4Z))*(z-M3Z)+M3Y)){sok=1; s++;}
				
			if(sok==1 && rok==1){     //反射かつスリットok
				while(z<Z0+BETA+GAMMA+NU){
					time_evolution(&t, &x, &y, &z ,vx ,vy ,vz);
				}
			}
				
			if(sok==1 && rok==0){     //反射しないかつスリットok
				while(z<Z0+BETA+GAMMA+NU){
					time_evolution(&t, &x, &y, &z ,vx ,vy ,vz);
				}
			}
		}
//粒子をテキストファイルへ		
			
		fprintf(f,"%lf %lf %lf %d %d %d\n",x,y,z,msok,sok,rok);
		
		
	}
	fprintf(f,"#beamnumber=%d\n",beamnumber);
    fprintf(f,"#time=%d\n",time);
	fprintf(f,"#ms=%d\n",ms);
	fprintf(f,"#s=%d\n",s);
	fprintf(f,"#ref=%d\n",ref);
	
	
	
	
	
		
	cout << "beamnumber=" << beamnumber << endl;
	cout << "time=" << time << endl;
	cout << "ms=" << ms << endl;
	cout << "s=" << s << endl;
	cout << "ref=" << ref << endl;
	
	
	fclose(f);
}
	
