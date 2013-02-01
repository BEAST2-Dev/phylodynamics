package beast.phylodynamics.util;

public class SEIRStateDouble {
	
    public double S;
    public double E;
    public double I;
    public double R;
    public double time;
    
    public SEIRStateDouble(double S, double E, double I, double R, double time) {

        this.S = S;
        this.E = E;
        this.I = I;
        this.R = R;
        this.time = time;
        
    }
    
    public SEIRStateDouble copy() {
    	return new SEIRStateDouble(S, E, I, R, time);
    }

}
