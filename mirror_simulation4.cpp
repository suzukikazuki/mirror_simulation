#include <iostream>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <math.h>

using namespace std;

const int N=300;   // kuans_energy�̕����̐�(����ɕύX�s��)
const int TRY=100000; // ���s��
const double PI=3.141592;
const double DELTA=5.0; //�~���[�X���b�g�̒Z�ӂ̒���(mm)
const double IPSI=60.0; //�~���[�X���b�g�̒��҂̒���(mm)
const double BETA=200; //�Օ��̂̒���(mm)
const double GAMMA=70.0; //�~���[�̒���(mm)
const double Z0=-(BETA+GAMMA); //z����0�_�̈ʒu(mm)
const double NU=2000; //�~���[���o����̒���(mm)
const int PERBEAM=10; //���ۂ̒����q���̂P�b������̗��q��
const int MESTIME=100000; //���莞��(s)
const double MKG=1.67*1e-27; //�����q�̎���(kg)
const double MMEV=939.6; //�����q�̎���(Mev)
const double TOFM=1240; //TOF���̒���(mm)
const double DT=0.0000001; //���Ԕ��W�̍��ݕ�
const double XMIN=-0.5*IPSI; //�X���b�g���ɂ������(�ǂɓ�����Ȃ�����)
const double XMAX=0.5*IPSI;
const double YMIN=-0.5*DELTA;
const double YMAX=0.5*DELTA;
const double YZMAX=PI/30; //0.0128
const double XZMAX=PI/30; //0.0205
const double DALPHA=PI/180; //�~���[�̌X��(z�����炙���ւ̊p)
const double M1Y=YMAX;  //�X������̃~���[��yz���W
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
	ifstream ifs("kuans_energy.txt"); // kuans_energy�̓ǂݍ���
	
	srand((unsigned int)time(NULL));
	
	
	
	
	int k[N];         //�K�v�Ȃ������t�@�C��kuans_energy.txt�̍����������̂��ʓ|�Ȃ̂ł���ɓ���Ă���
	int a;       //E�������̗���(0~299)
	double b[N];      //kuans_energy�̕��z
	int tof;     //����i�ł�TOF(0~3000mus)
	int d;       //counts���̗���(0~399)
	double azi;    //�Ɋp
	double pol;    //���ʊp
	double x, y, z;  //���W
	double vx,vy,vz; //���x����
	int beamnumber=0; //���q��
	int time=0; //����
	double v; //�����q�̑��x(mm/s)
	double lambda;  //�����q�̔g��(�I���O�X�g���[��)
	double ene; //�����q�̃G�l���M�[(MeV)
	double t=0; //���Ԕ��W�ɂ�����������
	int ms=0; //�~���[�X���b�g�ɓ�������
	int ref=0; //���˂��ꂽ�����q�̐�
	int s=0; //�X���b�g�𔲂�����
	double sine; //�p�x�����������̏d�݂Â��ɕK�v�ȕϐ�
	double vh; //���x�̐�������(YZ����)
	double dbeta; //�ŏ��̗��q�̓��˕���(z�����炙���ւ̊p)
	
	for(int i=0; i<N; i++){
		ifs >> k[i] >> b[i];
	}
	
	for(int i=0; i<TRY; i++){
		
		int msok=0; //�~���[�X���b�g�ɓ��������ǂ���
		int rok=0; //���˂������ǂ���
		int sok=0; //�~���[�𔲂�����
		
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

// �Օ��̂̓���
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
			
// �~���[�̓���			
		if(msok==1){
			dbeta=atan(vy/vz);
			vh=sqrt(v*v-vx*vx);
			while(XMIN<x && x<XMAX && y<(((M1Y-M4Y)/(M1Z-M4Z))*(z-M1Z)+M1Y) && y<(((M3Y-M4Y)/(M3Z-M4Z))*(z-M3Z)+M3Y) && y>(((M2Y-M3Y)/(M2Z-M3Z))*(z-M2Z)+M2Y)){
					
				time_evolution(&t, &x, &y, &z ,vx ,vy ,vz);
				if(y<(((M2Y-M3Y)/(M2Z-M3Z))*(z-M2Z)+M2Y) && lambda>3 && 0<DALPHA-dbeta && DALPHA-dbeta<PI/180){
					rok=1; ref++; goto EXIT;  //����
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
				
			if(sok==1 && rok==1){     //���˂��X���b�gok
				while(z<Z0+BETA+GAMMA+NU){
					time_evolution(&t, &x, &y, &z ,vx ,vy ,vz);
				}
			}
				
			if(sok==1 && rok==0){     //���˂��Ȃ����X���b�gok
				while(z<Z0+BETA+GAMMA+NU){
					time_evolution(&t, &x, &y, &z ,vx ,vy ,vz);
				}
			}
		}
//���q���e�L�X�g�t�@�C����		
			
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
	
