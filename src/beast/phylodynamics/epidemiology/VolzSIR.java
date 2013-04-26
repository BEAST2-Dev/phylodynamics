package beast.phylodynamics.epidemiology;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.coalescent.PopulationFunction;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Alexei Drummond
 */
@Description("Population function based on Volz 1999 `coalescent' likelihood.")
public class VolzSIR extends PopulationFunction.Abstract implements Loggable {
    public Input<RealParameter> n_S_Parameter = new Input<RealParameter>("n_S0",
            "the number of susceptibles at time of origin (defaults to 1000). ");
    public Input<RealParameter> betaParameter = new Input<RealParameter>("beta",
            "the mass action rate of infection.");
    public Input<RealParameter> gammaParameter = new Input<RealParameter>("gamma",
            "the per-infected rate of recovery.");
    public Input<RealParameter> originParameter = new Input<RealParameter>("origin",
            "the time before the root that the first infection occurred.");

    public Input<Double> integrationStepInput = new Input<Double>("integrationStep",
            "length of integration time step.", Validate.REQUIRED);
    public Input<Double> finishingThresholdInput = new Input<Double>("finishingThreshold",
            "Integration will finish when infected pop drops below this.", 1.0);
    
    public Input<Integer> statesToLogInput = new Input<Integer>("statesToLog",
            "Number of states to log. (Default 100.)", 100);
    
    //
    // Private stuff
    //
    
    private boolean dirty;
    private List<Double> effectivePopSize;
    private List<Double> NStraj, NItraj;
    
    //
    // Public stuff
    //
    
    @Override
    public void initAndValidate() throws Exception {
        if (betaParameter.get() != null) {
            betaParameter.get().setBounds(
            		Math.max(0.0, betaParameter.get().getLower()),
                    betaParameter.get().getUpper());
        }
        if (gammaParameter.get() != null) {
            gammaParameter.get().setBounds(
                    Math.max(0.0, gammaParameter.get().getLower()),
                    gammaParameter.get().getUpper());
        }
        if (originParameter.get() != null) {
            originParameter.get().setBounds(
                    Math.max(0.0, originParameter.get().getLower()),
                    originParameter.get().getUpper());
        }
        if (n_S_Parameter.get() != null) {
            n_S_Parameter.get().setBounds(
                    Math.max(1.0, n_S_Parameter.get().getLower()),
                    n_S_Parameter.get().getUpper());
        }

        effectivePopSize = new ArrayList<Double>();
        NStraj = new ArrayList<Double>();
        NItraj = new ArrayList<Double>();
        
        dirty = true;
    }
    
    /**
     * Update deterministic trajectory.
     */
    private void update() {
        if (!dirty)
            return;
        
        // Short-hand for model parameters:
        double beta = betaParameter.get().getValue();
        double gamma = gammaParameter.get().getValue();
        double NS0 = n_S_Parameter.get().getValue();
        
        double dt = integrationStepInput.get();
        double threshNI = finishingThresholdInput.get();
        
        // Clear old trajectory
        effectivePopSize.clear();
        NStraj.clear();
        NItraj.clear();
        
        // Set up initial conditions:
        double NS = NS0;
        double NI = 1.0;

        effectivePopSize.add(NI/(2.0*beta*NS));
        NStraj.add(NS);
        NItraj.add(NI);

        // Evaluate derivatives:
        double dNSdt, dNIdt = 0.0;
        do {
            
            double NSmid = NS;
            double NImid = NI;
            for (int iter=0; iter<3; iter++) {
                dNSdt = -beta*NSmid*NImid;
                dNIdt = beta*NSmid*NImid - gamma*NImid;
                
                NSmid = NS + 0.5*dt*dNSdt;
                NImid = NI + 0.5*dt*dNIdt;
            }
            NS = 2.0*NSmid - NS;
            NI = 2.0*NImid - NI;
            
            effectivePopSize.add(NI/(2.0*beta*NS));
            NStraj.add(NS);
            NItraj.add(NI);
        } while (dNIdt>0 || NI>threshNI);
        
        dirty = false;
    }
    

    // Implementation of abstract methods

    @Override
    public List<String> getParameterIds() {

        String[] parameterIds = new String[] {
            betaParameter.get().getID(),
            originParameter.get().getID()

        };
        return Arrays.asList(parameterIds);

    }

    /**
     * "Effective" population size for calculating coalescent likelihood
     * under Volz's model.
     * 
     * @param t
     * @return NI(t)/(2*beta*NS(t))
     */
    @Override
    public double getPopSize(double t) {
        update();
        
        // Choose which index into integration lattice to use:        
        double tForward = originParameter.get().getValue() - t;
        int tidx = (int)Math.floor(tForward/integrationStepInput.get());

        // Use initial or final state of trajectory if t outside the bounds of the
        // simulation.  This is a CLUDEGE to deal with trees which don't fit
        // the trajectories at all.
        if (tidx>=effectivePopSize.size())
            tidx = effectivePopSize.size()-1;
        else if (tidx<0)
            tidx = 0;

        return effectivePopSize.get(tidx);
    }
    
        
    @Override
    public double getIntegral(double start, double finish) {
        return getNumericalIntegral(start, finish);
    }

    @Override
    public double getIntensity(double t) {
        throw new UnsupportedOperationException("Not supported yet.");
        // only needed if not using numerical integration
    }

    @Override
    public double getInverseIntensity(double x) {
        throw new UnsupportedOperationException("Not supported yet.");
        // only needed for simulation
    }

    /*
     * CalculationNode interface
     */
    
    @Override
    public boolean requiresRecalculation() {
        dirty = true;
        return true;
    }
    
    @Override
    public void restore() {
        dirty = true;
        super.restore();
    }
    
    /*
     * Loggable interface
     */

    @Override
    public void init(PrintStream out) throws Exception {
        for (int i=0; i<statesToLogInput.get(); i++) {
            out.format("S%d\t",i);
            out.format("I%d\t",i);
            out.format("t%d\t",i);
        }
    }

    @Override
    public void log(int nSample, PrintStream out) {
        double dt = integrationStepInput.get();
        double tend = NStraj.size()*dt;
        
        double delta = tend/(statesToLogInput.get()-1);
        
        for (int i=0; i<statesToLogInput.get(); i++) {
            
            double t = delta*i;
            int tidx = (int)Math.round(t/dt);
            
            if (tidx>=NStraj.size())
                tidx = NStraj.size()-1;
            
            out.format("%g\t", NStraj.get(tidx));
            out.format("%g\t", NItraj.get(tidx));
            out.format("%g\t", t);
        }
    }

    @Override
    public void close(PrintStream out) {
    }
}