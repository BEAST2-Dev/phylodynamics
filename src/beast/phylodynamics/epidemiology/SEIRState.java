package beast.phylodynamics.epidemiology;

public class SEIRState {

    public double S;
    public double E;
    public double I;
    public double R;
    public double time;

    public SEIRState(double S, double E, double I, double R, double time) {

        this.S = S;
        this.E = E;
        this.I = I;
        this.R = R;
        this.time = time;
    }
    
    public SEIRState copy() {
    	return new SEIRState(S, E, I, R, time);
    }
}
