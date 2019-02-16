
import java.util.Arrays;
import edu.emory.mathcs.jtransforms.fft.*;


//import edu.emory.mathcs.jtransforms.fft.DoubleFFT_3D;

public class ComplexTensor3 {
	double[][][] el;
	int I,J,K;
	
	ComplexTensor3(){}
	ComplexTensor3(int I,int J,int K){
		el=new double[I][J][2*K];
		this.I=I;
		this.J=J;
		this.K=K;
	}

	public ComplexTensor3 add(ComplexTensor3 A){
		if(I!=A.I || J!=A.J || K!=A.K) throw new IllegalArgumentException("Tensor dimensions do not agree");
		ComplexTensor3 R=new ComplexTensor3(I,J,K);
		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				for(int k=0;k<2*K;k++)
					R.el[i][j][k]=el[i][j][k]+A.el[i][j][k];
		return R;
	}

	public ComplexTensor3 sub(ComplexTensor3 A){
		if(I!=A.I || J!=A.J || K!=A.K) throw new IllegalArgumentException("Tensor dimensions do not agree");
		ComplexTensor3 R=new ComplexTensor3(I,J,K);
		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				for(int k=0;k<2*K;k++)
					R.el[i][j][k]-=A.el[i][j][k];
		return R;
	}
	
	public void fillWith(ComplexTensor3 A){
		if(I!=A.I || J!=A.J || K!=A.K) throw new IllegalArgumentException("Tensor dimensions do not agree");
		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				el[i][j]=Arrays.copyOf(A.el[i][j], 2*K);
			}
	
	public ComplexTensor3 dot121(ComplexTensor3 A){
		if(I!=A.I || J!=A.J || K!=A.K) throw new IllegalArgumentException("Tensor dimensions do not agree");
		ComplexTensor3 R=new ComplexTensor3(I,J,K);
		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				for(int k=0;k<K;k++){
					R.el[i][j][2*k]=el[i][j][2*k]*A.el[i][j][2*k]-el[i][j][2*k+1]*A.el[i][j][2*k+1];
					R.el[i][j][2*k+1]=el[i][j][2*k]*A.el[i][j][2*k+1]+el[i][j][2*k+1]*A.el[i][j][2*k];
				}
		return R;
			}
	
	public void times(double a){
	
		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				for(int k=0;k<2*K;k++)
					el[i][j][k]*=a;
	}
	
	public void ones(){
		
		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				for(int k=0;k<K;k++){
					el[i][j][2*k]=1;
					el[i][j][2*k+1]=0;
									}
		}
	
	public void FFT(){
		DoubleFFT_3D ft3=new DoubleFFT_3D(I,J,K);
		ft3.complexForward(el);
	}
	
	public void iFFT(){
		DoubleFFT_3D ft3=new DoubleFFT_3D(I,J,K);
		ft3.complexInverse(el,true);
	}

}
