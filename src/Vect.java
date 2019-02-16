
import static java.lang.Math.*;

import java.util.Arrays;
public class Vect {
	Vect(){};
	double[] el;
	int length;
	
	public Vect(int I){
		this.length=I;
		this.el=new double[I];
	}
	
	public void set(double[] u){
		this.el=u;
		this.length=u.length;		
	}

	public void set(int[] u){
		for(int i=0;i<u.length;i++)
		this.el[i]=(double)u[i];
		this.length=u.length;
	}
	
	public Vect add(Vect v){
	
		if(this.length!=v.length) throw new NullPointerException("vectrs have different lengths");
		Vect w=new Vect(v.length);
		for(int i=0;i<v.length;i++)
			w.el[i]=this.el[i]+v.el[i];
		return w;
	}
	
	public Vect sub(Vect v){
		
		if(this.length!=v.length) throw new NullPointerException("vectrs have different lengths");
		Vect w=new Vect(v.length);
		for(int i=0;i<v.length;i++)
			w.el[i]=this.el[i]-v.el[i];
		return w;
	}
	
	public Vect rand(int I){
		
		Vect v=new Vect(I);
		for(int i=0;i<I;i++)
			v.el[i]=random();
		return v;
	}
	public Vect rand(){
		
		Vect v=new Vect(length);
		for(int i=0;i<length;i++)
			v.el[i]=random();
		return v;
	}
	
	public Vect ones(int I){
		
		Vect v=new Vect(I);
		for(int i=0;i<I;i++)
			v.el[i]=1;
		return v;
	}

	public Vect flr(){
		int I=length;
		Vect v=new Vect(I);
		for(int i=0;i<I;i++)
			v.el[i]=floor(el[i]);
		return v;
	}
	public Vect rand(int I,double a, double b){
		
		Vect v=new Vect(I);
		for(int i=0;i<I;i++)
			v.el[i]=a+(b-a)*random();
		return v;
	}
	
public Vect rand(int I,int a, int b){
		
		Vect v=new Vect(I);
		for(int i=0;i<I;i++)
			v.el[i]=a+(b-a)*random();
		return v;
	}


	public Vect linspace(double a, double b,int I){
		
		Vect v=new Vect(I);
		double d=(b-a)/(I-1);
		for(int i=0;i<I;i++)
			v.el[i]=a+i*d;
		return v;
	}
	
	public Vect sqspace(double a, double b,int I){
		
		Vect v=new Vect(I);
		double d=(b-a)/(I-1);
		v.el[0]=0;
		for(int i=1;i<I;i++)
			v.el[i]=v.el[i-1]+1+0.2*abs(i-(double)I/2);
		v=v.time((b-a)/v.el[I-1]);
		for(int i=0;i<I;i++)
			v.el[i]+=a;
		return v;
	}
	
public Vect randspace(double a, double b,int I,double r){
		
		Vect v=new Vect();
		double d=(b-a)/(I-1);
		v=v.rand(I,-r*d,r*d);	
		v=v.add(v.linspace(a,b,I));
		v.el[0]=a;
		v.el[I-1]=b;

		return v;
	}
	
	public Vect time(double a){
		
		Vect v=new Vect(this.length);
		for(int i=0;i<this.length;i++)
			v.el[i]=a*this.el[i];
		return v;
	}
	public Vect times(double a){
		
		Vect v=new Vect(this.length);
		for(int i=0;i<this.length;i++)
			v.el[i]=a*this.el[i];
		return v;
	}

	public Vect times(int a){
		
		Vect v=new Vect(this.length);
		for(int i=0;i<this.length;i++)
			v.el[i]=a*this.el[i];
		return v;
	}

	
	public double dot(Vect u){
		
		if(this.length!=u.length) throw new NullPointerException("vectrs have different lengths");
		double s=0;
		for(int i=0;i<u.length;i++)
			s=s+this.el[i]*u.el[i];
		return s;
	}
	
	public Vect cross(Vect u){
		if(length!=u.length) throw new NullPointerException("vectrs have different lengths");
		if(length>3) throw new NullPointerException("Cross product is not defined for dimentions higher than 3");
		Vect w=new Vect(3);
			w.el[0]=el[1]*u.el[2]-el[2]*u.el[1];
			w.el[1]=el[2]*u.el[0]-el[0]*u.el[2];
			w.el[2]=el[0]*u.el[1]-el[1]*u.el[0];
			return w;
		}

	public double norm(){
		double s=0;
		for(int i=0;i<this.length;i++)
			s=s+this.el[i]*this.el[i];
		return sqrt(s);
	}
	
	public  double[] angles(Vect u){
		if(length>3) throw new NullPointerException("Vector dimentions is higher than 3");

		double[] ThetPhi=new double[2];

		ThetPhi[0]=acos(u.el[2]);
		if(u.el[0]==0) 
			ThetPhi[1]=(2-signum(u.el[1]))*PI/2;
		else if(u.el[1]==0) 
			ThetPhi[1]=(1-signum(u.el[0]))*PI/2;
		else if(u.el[0]*u.el[1]<0) 
			ThetPhi[1]=(3-signum(u.el[1]))*PI/2+atan(u.el[1]/u.el[0]);
		else if(u.el[0]*u.el[1]>0) 
			ThetPhi[1]=(1-signum(u.el[0]))*PI/2+atan(u.el[1]/u.el[0]);
		
		return ThetPhi;
	}

	
}
