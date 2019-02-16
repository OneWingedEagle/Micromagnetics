
import static java.lang.Math.*;
import java.util.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.*;

import javax.swing.BorderFactory;
import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;

public class Model {
	private ComplexTensor3 DTXX,DTXY,DTXZ,DTYY,DTYZ,DTZZ;
	public int I,J,K;
	public double Lx,Ly,Lz,dxnm,dynm,dznm,dx,dy,dz;
	public int Nx,Ny,Nz;
	private int dim=3;
	private double u0=4*PI*1e-7;
	public double mu0Ms;
	private double drx,dry,drz,dv;
	public double A,Ms,K1,K2,Eexch=0,Eext=0,Eani=0,Ed=0,Etotal=0;
	public double alpha,gamma,betta,dt;
	public  double AngMax;
	private  Cell[][][] cell;
	public Vect uniAniAxis=new Vect(3);
	public boolean Rot,LLG,PreviousInitial,ModifEex,ModifEd,CubicAni,done;




Model(){}

	public void setCells(){
		dx=dxnm*1e-9;
		dy=dynm*1e-9;
		dz=dznm*1e-9;
		dv=dx*dy*dz;
		cell=new Cell[I][J][K];
		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				for(int k=0;k<K;k++)
					cell[i][j][k]=new Cell();
	}
	
	//==================================== This method  initializes M from the previous run.

	public void setInitialPrev(String prevMagFile){

		int i,j,k;
		try{
		FileReader fr = new FileReader(prevMagFile);
		Scanner scr = new Scanner(fr);
		Lx=scr.nextDouble();
		Ly=scr.nextDouble();
		Lz=scr.nextDouble();
		I=scr.nextInt();
		J=scr.nextInt();
		K=scr.nextInt();
		for(i=0;i<I;i++)
			for(j=0;j<J;j++)
				for(k=0;k<K;k++){
					cell[i][j][k].m.el[0]=scr.nextDouble();
					cell[i][j][k].m.el[1]=scr.nextDouble();
					cell[i][j][k].m.el[2]=scr.nextDouble();	

				}
		scr.close();
		fr.close();
		}
		catch(IOException e){System.err.println("Error in reading m from the previous run");
		}
	
	}
	
	//===================================== This method  initializes M vith new valuse.
	public void setInitialNew(double Theta, double ThetaRange,double Phi, double PhiRange){


		double thetaMean=Theta/180*PI;
		double thetaRange=ThetaRange/180*PI;
		double phiMean=Phi/180*PI;
		double phiRange=PhiRange/180*PI;
		double theta,phi;
		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				for(int k=0;k<K;k++){
					theta=thetaMean+thetaRange*(.5-random());
					phi=phiMean+phiRange*(.5-random());
					cell[i][j][k].m.el[0]=sin(theta)*cos(phi);
					cell[i][j][k].m.el[1]=sin(theta)*sin(phi);
					cell[i][j][k].m.el[2]=cos(theta);
					cell[i][j][k].m=cell[i][j][k].m.time(1.0/cell[i][j][k].m.norm());
					
					

				}	
		}
		
	

	//====================================================This method calculates  FFT of m components.
	public ComplexTensor3 getFFTofmComponent(int p){
		ComplexTensor3 mp=new ComplexTensor3(2*I,2*J,2*K);
		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				for(int k=0;k<K;k++)
					mp.el[i][j][2*k]=cell[i][j][k].m.el[p];
		mp.FFT();

		return mp;	
	}

	//=====================================================This method calculates  the demagnetizing  energy and field.

	public void calHd(){
		ComplexTensor3 MX=getFFTofmComponent(0);
		ComplexTensor3 MY=getFFTofmComponent(1);
		ComplexTensor3 MZ=getFFTofmComponent(2);
		ComplexTensor3 HU=new ComplexTensor3();

		HU=DTXX.dot121(MX).add(DTXY.dot121(MY)).add(DTXZ.dot121(MZ));
		HU.iFFT();
		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				for(int k=0;k<K;k++)
					cell[i][j][k].Hd.el[0]=-Ms*HU.el[i+I][j+J][2*k+2*K]/(32*PI);

		HU=DTXY.dot121(MX).add(DTYY.dot121(MY)).add(DTYZ.dot121(MZ));
		HU.iFFT();

		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				for(int k=0;k<K;k++)
					cell[i][j][k].Hd.el[1]=-Ms*HU.el[i+I][j+J][2*k+2*K]/(32*PI);

		HU=DTXZ.dot121(MX).add(DTYZ.dot121(MY)).add(DTZZ.dot121(MZ));
		HU.iFFT();

		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				for(int k=0;k<K;k++)
					cell[i][j][k].Hd.el[2]=-Ms*HU.el[i+I][j+J][2*k+2*K]/(32*PI);
		Ed=0;
		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				for(int k=0;k<K;k++)
					Ed+=cell[i][j][k].m.dot(cell[i][j][k].Hd);
		Ed*=-Ms/2*u0*dv;

	}
	
	//======================================================== This method calculates the exchange energy and field.

	public void calHexch()
	{
		Eexch=0;
		double Ax=2*A/mu0Ms;
		Vect[] NN;
		Vect SN=new Vect(dim);
		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				for(int k=0;k<K;k++){
					NN=NN6(i,j,k);
					Vect mcell=cell[i][j][k].m;
					for(int p=0;p<6;p++)
						SN=SN.add(NN[p]);
					Eexch+=(2-mcell.dot(NN[0])-mcell.dot(NN[1]))/(dx*dx)+(2-mcell.dot(NN[2])-mcell.dot(NN[3]))/(dy*dy)+(2-mcell.dot(NN[4])-mcell.dot(NN[5]))/(dz*dz);
					cell[i][j][k].Hexch=NN[0].add(NN[1]).times(1.0/(dx*dx)).add(NN[2].add(NN[3]).times(1.0/(dy*dy))).add(NN[4].add(NN[5]).times(1.0/(dz*dz))).times(Ax);;
				}
		Eexch*=A*dv;
	}

	//========================================================This method calculates the anisotropy energy and field.

	public void calHani(){

		Vect mcell=new Vect(dim);
		double mdn=0;
		
		Eani=0;
		if(CubicAni)
		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				for(int k=0;k<K;k++){
					mcell=cell[i][j][k].m;
					Eani=Eani+K1*(pow(mcell.el[0]*mcell.el[1],2)+pow(mcell.el[0]*mcell.el[2],2)+pow(mcell.el[1]*mcell.el[2],2))+K2*pow(mcell.el[0]*mcell.el[1]*mcell.el[2],2);
					cell[i][j][k].Hani.el[0]=-2.0*K1*mcell.el[0]*(1-pow(mcell.el[0],2))-2.0*K2*pow(mcell.el[1]*mcell.el[2],2)*mcell.el[0]/mu0Ms;
					cell[i][j][k].Hani.el[1]=-2.0*K1*mcell.el[1]*(1-pow(mcell.el[1],2))-2.0*K2*pow(mcell.el[0]*mcell.el[2],2)*mcell.el[1]/mu0Ms;
					cell[i][j][k].Hani.el[2]=-2.0*K1*mcell.el[2]*(1-pow(mcell.el[2],2))-2.0*K2*pow(mcell.el[0]*mcell.el[1],2)*mcell.el[2]/mu0Ms;

				}
		else
			for(int i=0;i<I;i++)
				for(int j=0;j<J;j++)
					for(int k=0;k<K;k++){
						mcell=cell[i][j][k].m;
						mdn=mcell.dot(uniAniAxis);
						Eani=Eani+K1*(1-pow(mdn,2))+K2*pow(1-pow(mdn,2),2);
						cell[i][j][k].Hani=uniAniAxis.times(2*mdn*(K1+2*K2*(1-pow(mdn,2)))/mu0Ms);

					}		
		Eani*=dv;
	}
	
	
	//===================================================This method sets the given applied field to all the cells.

	public void setHext(double Ha,double thetaHa,double phiHa){
		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				for(int k=0;k<K;k++){
					cell[i][j][k].Hext.el[0]=Ha*sin(thetaHa)*cos(phiHa);
					cell[i][j][k].Hext.el[1]=Ha*sin(thetaHa)*sin(phiHa);
					cell[i][j][k].Hext.el[2]=Ha*cos(thetaHa);
					
				}

	}

	//==================================================== This method calculates Zeeman energy.

	public void calEext(){
		Eext=0;
		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				for(int k=0;k<K;k++)
					Eext+=cell[i][j][k].m.dot(cell[i][j][k].Hext);
		Eext*=-Ms*u0*dv;
	}

	//================================================== This method adds up the energies to get the total energy.

	public void calEtotal(){
		Etotal=Eexch+Eext+Ed+Eani;		
	}
	
	public void movem(){
		double dt1=dt*Ms*gamma;
		double gammaL=1/(1+pow(alpha,2));
		double lambda=alpha*gammaL;
		AngMax=0;
		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				for(int k=0;k<K;k++){
					 cell[i][j][k].calAngDiff();
	                    if(cell[i][j][k].angdiff > AngMax)
	                    	  AngMax = cell[i][j][k].angdiff;
					if(LLG)
					cell[i][j][k].LLG(gammaL, lambda, Ms, dt1);
					else
					cell[i][j][k].rotation(betta);
					
					}
	}

	//=============================================== This method gets the 6 nearest neighbor required in calculation of the exchange energy.

	public  Vect[] NN6(int i,int j, int k){

		Vect[] NN=new Vect[6];
		for(int p=0;p<6;p++)
			NN[p]=new Vect(dim);

		if(j==0){ for(int p=0;p<3;p++) NN[0].el[p]=cell[i][j][k].m.el[p];}
		else NN[0]=cell[i][j-1][k].m;

		if(j==J-1){ for(int p=0;p<3;p++) NN[1].el[p]=cell[i][j][k].m.el[p];}
		else NN[1]=cell[i][j+1][k].m;

		if(i==0){ for(int p=0;p<3;p++) NN[2].el[p]=cell[i][j][k].m.el[p];}
		else NN[2]=cell[i-1][j][k].m;

		if(i==I-1){ for(int p=0;p<3;p++) NN[3].el[p]=cell[i][j][k].m.el[p];}
		else NN[3]=cell[i+1][j][k].m;

		if(k==0){ for(int p=0;p<3;p++) NN[4].el[p]=cell[i][j][k].m.el[p];}
		else NN[4]=cell[i][j][k-1].m;

		if(k==K-1) { for(int p=0;p<3;p++) NN[5].el[p]=cell[i][j][k].m.el[p];}
		else NN[5]=cell[i][j][k+1].m;
		return NN;
	}

	//================================================== This method calculates the effective field by summing the previously calculated fields.

	public void calHeff(){

		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				for(int k=0;k<K;k++)			
					cell[i][j][k].Heff=cell[i][j][k].Hd.add(cell[i][j][k].Hexch).add(cell[i][j][k].Hext).add(cell[i][j][k].Hani);
	}
	

	
	//===================================================== This calculates FFT of the demagnetizing tensor.
	public  void setDemagTensorFFT(){
		int I2,J2,K2;
		I2=2*I;J2=2*J;K2=2*K;
		int counter=1;
		int total=I2*J2*K2;
		drx=1; dry=dy/dx; drz=dz/dx;
	
	    System.out.println(" Calculating demagnetizing tensors ... ");
	    System.out.println();
	    System.out.println(" Progress... ");
	    System.out.println();

		DTXX=new ComplexTensor3(I2,J2,K2);
		DTXY=new ComplexTensor3(I2,J2,K2);
		DTXZ=new ComplexTensor3(I2,J2,K2);
		DTYY=new ComplexTensor3(I2,J2,K2);
		DTYZ=new ComplexTensor3(I2,J2,K2);
		DTZZ=new ComplexTensor3(I2,J2,K2);
		double v=drx*dry*drz;
		for(int q=0;q<I2;q++)
			for(int p=0;p<J2;p++)
				for(int r=0;r<K2;r++){
					counter++;
					if(counter%(total/10)==0) {
						System.out.println(" ... "+10*counter/(total/10)+" %");
					}

					for(int i=1;i<4;i++)
						for(int j=1;j<4;j++)
							for(int k=1;k<4;k++){
								DTXX.el[q][p][2*r]+=pow(-1,i+j+k-1)*sn(i)*sn(j)*sn(k)*F1(drx*(p-J)+a(j,drx),dry*(q-I)+a(i,dry),drz*(r-K)+a(k,drz));
								DTYY.el[q][p][2*r]+=pow(-1,i+j+k-1)*sn(i)*sn(j)*sn(k)*F1(dry*(q-I)+a(i,dry),drz*(r-K)+a(k,drz),drx*(p-J)+a(j,drx));
								DTZZ.el[q][p][2*r]+=pow(-1,i+j+k-1)*sn(i)*sn(j)*sn(k)*F1(drz*(r-K)+a(k,drz),drx*(p-J)+a(j,drx),dry*(q-I)+a(i,dry));
								DTXY.el[q][p][2*r]+=pow(-1,i+j+k-1)*sn(i)*sn(j)*sn(k)*F2(drz*(r-K)+a(k,drz),drx*(p-J)+a(j,drx),dry*(q-I)+a(i,dry));
								DTXZ.el[q][p][2*r]+=pow(-1,i+j+k-1)*sn(i)*sn(j)*sn(k)*F2(dry*(q-I)+a(i,dry),drx*(p-J)+a(j,drx),drz*(r-K)+a(k,drz));
								DTYZ.el[q][p][2*r]+=pow(-1,i+j+k-1)*sn(i)*sn(j)*sn(k)*F2(drx*(p-J)+a(j,drx),dry*(q-I)+a(i,dry),drz*(r-K)+a(k,drz));
							}

					DTXX.el[q][p][2*r]*=8.0/v;
					DTYY.el[q][p][2*r]*=8.0/v;
					DTZZ.el[q][p][2*r]*=8.0/v;
					DTXY.el[q][p][2*r]*=8.0/v;
					DTXZ.el[q][p][2*r]*=8.0/v;
					DTYZ.el[q][p][2*r]*=8.0/v;

				}
		counter=1;

		DTXX.FFT();	
		DTXY.FFT();
		DTXZ.FFT();
		DTYY.FFT();
		DTYZ.FFT();
		DTZZ.FFT();
		System.out.println("Calculation of the demagnetizing factor completed.");

	}

	//===================== function required in the calculation of the demagnetizing field;
	public static double a(int i, double drx){
		if(i==1) return -drx;
		else if(i==2) return 0;
		else return drx;
	}

	//===================== function required in the calculation of the demagnetizing field;
	public static int sn(int i){
		int sn=1;

		if(i==2) sn=2;

		return sn;
	}


	//===================== function required in the calculation of the demagnetizing field;
	public static double F1(double X, double Y, double Z){
		double R,F1;

		R=sqrt(pow(X,2)+pow(Y,2)+pow(Z,2));

		if(R==0)
			F1=0;
		else if(X==0 && Y==0 && Z!=0)
			F1=Z*Z*abs(Z)/6;


		else if(X==0 && Y!=0 && Z!=0)
			F1=.5*Y*pow(Z,2)*log(R-Y)
			+.5*Z*pow(Y,2)*log(R-Z)
			+(Y*Y+Z*Z)*R/6;

		else if(X==0 && Y!=0 && Z==0)
			F1=-.5*Y*pow(X,2)*log(R)+(Y*Y-2*X*X)*R/6;

		else
			F1=X*Y*Z*atan(Y*Z/(X*R))+
			.5*Y*(pow(Z,2)-pow(X,2))*log(R-Y)
			+.5*Z*(pow(Y,2)-pow(X,2))*log(R-Z)
			+(Y*Y+Z*Z-2*X*X)*R/6;

		return F1;
	}

	//===================== function required in the calculation of the demagnetizing field;
	public static double F2(double Z, double X, double Y){
		double R,F2;

		R=sqrt(pow(X,2)+pow(Y,2)+pow(Z,2));

		if(X!=0 && Y!=0 && Z!=0)	
			F2=-X*Y*Z*log(R+Z)+
			Y*(pow(Y,2)-3*pow(Z,2))*log(R+X)/6.0
			+X*(pow(X,2)-3*pow(Z,2))*log(R+Y)/6.0
			+.5*X*X*Z*atan(Y*Z/(X*R))
			+.5*Y*Y*Z*atan(X*Z/(Y*R))
			+pow(Z,3)/6*atan(X*Y/(Z*R))+X*Y*R/3;

		else if(X==0 && Y!=0 && Z!=0)
			F2=Y*(pow(Y,2)-3*pow(Z,2))*log(R)/6.0;

		else if(X!=0 && Y==0 && Z!=0)
			F2=X*(pow(X,2)-3*pow(Z,2))*log(R)/6.0;

		else if(X==0 && Y==0 && Z!=0)
			F2=0;

		else if(X!=0 && Y!=0 && Z==0)
			F2=pow(Y,3)*log(R+X)/6.0+pow(X,3)*log(R+Y)/6.0+X*Y*R/3;

		else if(X==0 && Y!=0 && Z==0)
			F2=pow(Y,3)*log(R)/6.0;

		else if(X!=0 && Y==0 && Z==0)
			F2=pow(X,3)*log(R)/6.0;
		else
			F2=0;
		return F2;

	}

	//====================================== This method writes normalized magnetization, m, to the given path

	public void writemag(String magFile){

		
		try{ PrintWriter pw = new  PrintWriter(new  FileWriter(magFile));

		int i,j,k,D=3;
		pw.println("flux");
		pw.println("3");
		pw.println(I*J*K);

		for(k=0;k<K;k++)
			for(i=0;i<I;i++)	
				for(j=0;j<J;j++)	
				{ 
					if(i==I-1 && j==J-1 && k==0)
						for(int p=0;p<D;p++)	
							pw.format("%10.6f",0.0);
					else
						for(int p=0;p<D;p++)	
							pw.format("%10.6f",cell[i][j][k].m.el[p]);
					pw.println();
				}

		pw.close();
		}

		catch (IOException e) { } 


	}
	
	//====================================== This method writes normalized magnetization, m, to the given path for the next initialization if required

	public void writemForIni(String prevMagFile){
		try{ PrintWriter pw = new  PrintWriter(new  FileWriter(prevMagFile));
		int D=3;
		pw.format("%.1f\n",Lx);
		pw.format("%.1f\n",Ly);
		pw.format("%.1f\n",Lz);
		pw.format("%6d\n",J);
		pw.format("%6d\n",I);
		pw.format("%6d\n",K);


		for(int i=0;i<I;i++)	
			for(int j=0;j<J;j++)
				for(int k=0;k<K;k++)

				{ 
					for(int p=0;p<D;p++)	
						pw.format("%10.6f",cell[i][j][k].m.el[p]);	
					pw.println();	
				}
		pw.close();
		}

		catch (IOException e) { } 

	}

	//============================================ This method writes the effective field to the given path

	public void writeHeff(String HeffFile){
		try{ PrintWriter pw = new  PrintWriter(new  FileWriter(HeffFile));

		int i,j,k,D=3;
		pw.println("flux");
		pw.println("3");
		pw.println(I*J*K);
		for(k=0;k<K;k++)
			for(i=0;i<I;i++)	
				for(j=0;j<J;j++)	
				{ 
					for(int p=0;p<D;p++)	
						pw.format("%15.2f",cell[i][j][k].Heff.el[p]);
					pw.println();
				}

		pw.close();
		}

		catch (IOException e) { } 


	}

	//========================================= This method writes the effective field of the surface cells only, to the given path
	public void writeHeffSurf(String HeffSurfFile){
		try{ PrintWriter pw = new  PrintWriter(new  FileWriter(HeffSurfFile));

		int i,j,k,D=3;
		pw.println("flux");
		pw.println("3");
		pw.println(I*J*K);
		for(k=0;k<K;k++)
			for(i=0;i<I;i++)	
				for(j=0;j<J;j++)	
				{ 
					for(int p=0;p<D;p++)	
						if(i==0 || j==0 || k==K-1)
					pw.format("%15.2f",cell[i][j][k].Heff.el[p]);	
						else 
							pw.format("%15.2f",0.0);		
				
					pw.println();
				}

		pw.close();
		}

		catch (IOException e) { } 
		
		done=true;


	}

	//================================================ This method writes the normalized magnetization of the surface cells only, to the given path
	public void writemagSurf(String magSurfFile){
		try{ PrintWriter pw = new  PrintWriter(new  FileWriter(magSurfFile));

		int i,j,k,D=3;
		pw.println("flux");
		pw.println("3");
		pw.println(I*J*K);
		for(k=0;k<K;k++)
			for(i=0;i<I;i++)	
				for(j=0;j<J;j++)	
				{ 
					for(int p=0;p<D;p++)	
						if(i==0 || j==0 || k==K-1)
					pw.format("%10.6f",cell[i][j][k].m.el[p]);	
						else 
							pw.format("%10.6f",0.0);		
				
					pw.println();
				}

		pw.close();
		}

		catch (IOException e) { } 


	}


	//================================================ This method writes mesh file to the given path

	public void writeMesh(String MeshFile){

		try { PrintWriter pw = new  PrintWriter(new  FileWriter(MeshFile));
		int Nn,Nx,Nxy,Ne,frx,frxy;
		Nn=(J+1)*(I+1)*(K+1);
		Nx=J+1;
		Nxy=(I+1)*(J+1);
		Ne=J*I*K;
		pw.println("hexahedron");
		pw.println("//Number_of_Node");
		pw.println(Nn);
		pw.println("//Number_of_Element");
		pw.println(Ne);
		pw.println("//Number_of_Region");
		pw.println("1");
		pw.println("//Factor");
		pw.println("1");

		frx=1;
		frxy=1;
		while(frxy*Nxy<Nn)	{

			pw.print((frxy*Nxy+frx)+",");
			pw.print((frxy*Nxy+1+frx)+",");
			pw.print((frxy*Nxy+Nx+1+frx)+",");
			pw.print((frxy*Nxy+Nx+frx)+",");
			pw.print(((frxy-1)*Nxy+frx)+",");
			pw.print(((frxy-1)*Nxy+1+frx)+",");
			pw.print(((frxy-1)*Nxy+Nx+1+frx)+",");
			pw.println(((frxy-1)*Nxy+Nx+frx)+",");


			if ((frx+1)%Nx==0 &&(frx+Nx+1)%Nxy>0 )
				frx=frx+2;
			else if ((frx+Nx+1)%Nxy==0 )
			{	
				frxy=frxy+1;
				frx=1;
			}
			else frx=frx+1;

		}
		double d=1.0/I;

		for(int k=0;k<K+1;k++)
			for(int j=0;j<I+1;j++)
				for(int i=0;i<J+1;i++)

				{
					pw.format("%15.5f,",i*drx);
					pw.format("%15.5f,",j*dry);
					pw.format("%15.5f,\n",k*drz);

				}
		pw.format("1,");
		pw.format("%10d,",Ne);
		pw.format("%10.5f,\n",J*d);
		pw.close();

		}
		catch (IOException e) { } 
	}

}

