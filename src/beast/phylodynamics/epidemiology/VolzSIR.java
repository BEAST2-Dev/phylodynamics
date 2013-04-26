package beast.phylodynamics.epidemiology;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.coalescent.PopulationFunction;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Alexei Drummond
 */
@Description("Population function based on Volz 1999 `coalescent' likelihood.")
public class VolzSIR extends PopulationFunction.Abstract {
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
    
    //
    // Private stuff
    //
    
    private boolean dirty;
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
        NStraj.clear();
        NItraj.clear();
        
        // Set up initial conditions:
        double NS = NS0;
        double NI = 1.0;

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
        
        double beta = betaParameter.get().getValue();
        double tForward = originParameter.get().getValue() - t;
        
        // Is this sensible?
        if (tForward<0.0)
            return 1.0/(2.0*beta*n_S_Parameter.get().getValue());
        
        // Choose which index into integration lattice to use:
        int tidx = (int)Math.floor(tForward/integrationStepInput.get());

        // Use last final state of trajectory if t outside the bounds of the
        // simulation.  This is a CLUDEGE to deal with trees which don't fit
        // the trajectories at all.
        if (tidx>=NItraj.size())
            tidx = NItraj.size()-1;

        return NItraj.get(tidx)/(2.0*beta*NStraj.get(tidx));
        
//        if (tidx<NItraj.size())
//            return NItraj.get(tidx)/(2.0*beta*NStraj.get(tidx));
//        else
//            return 0.0;
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
}