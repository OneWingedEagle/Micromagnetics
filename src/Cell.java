

import static java.lang.Math.*;

public class Cell {
	
	int comp=3;
	public Vect m=new Vect(comp);
	public Vect Heff=new Vect(comp);
	public Vect Hd=new Vect(comp);
	public Vect Hexch=new Vect(comp);
	public Vect Hext=new Vect(comp);
	public Vect Hani=new Vect(comp);
    public double angdiff;


	Cell(){}
	
	public void rotation(double beta){
		double nH=Heff.norm();
		double eps=1e-3;
		Vect dm=new Vect();
		dm=dm.rand(comp,-eps,eps);

		if(nH>0 && m.norm()>.10) 
			{
		Vect h=Heff.times(1.0/nH);
		double mdoth=m.dot(h);
		if(mdoth<-0.8){
			if(mdoth==-1) m=m.add(dm);
			h=h.add(m);
			h=h.times(1.0/h.norm());
			}
		
			if(mdoth>.995)
				m=h;
			else{
		m=m.times(1-beta).add(h.times(beta));	
		m=m.times(1.0/m.norm());
			}
		
			}

		}

    	public  void LLG(double gammaL, double lambda, double Ms, double dt){
    		
    		double eps=1e-4;
    		
    		Vect dm=new Vect();
    		dm=dm.rand(comp,-eps,eps);
    		Vect heff=Heff.times(1.0/Ms);
    		double nheff=heff.norm();
    		
		if(nheff>0){	
				Vect mxh=m.cross(heff);   	
    			if(m.dot(heff)/nheff==-1) m=m.add(dm);
       			m=m.add(((mxh.times(-gammaL)).add(m.cross(mxh).times(-lambda))).times(dt));
    			double normm=m.norm();
    			if(normm>1.00001 || normm<.99999)
    		    m=m.times(1.0/normm);
    				
   			}
		
		
		
    	}
    	   public void calAngDiff()
    	    {
    	        double nHeff = Heff.norm();
    	        if(nHeff == 0.0D)
    	        {
    	            angdiff = 0.0D;
    	        } else
    	        {
    	            double mh = m.dot(Heff) / nHeff;
    	            if(mh >= 1.0D)
    	                angdiff = 0.0D;
    	            else
    	            if(mh <= -1D)
    	                angdiff = PI;
    	            else
    	                angdiff = Math.acos(mh);
    	        }
    	    }
}