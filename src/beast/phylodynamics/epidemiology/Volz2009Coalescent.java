package beast.phylodynamics.epidemiology;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.TreeDistribution;
import beast.evolution.tree.coalescent.IntervalType;
import beast.evolution.tree.coalescent.TreeIntervals;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.ClassicalRungeKuttaIntegrator;

/**
 * @author Alexei Drummond
 */

@Description("Calculates the Volz-2009 probability of a beast.tree.")
public class Volz2009Coalescent extends TreeDistribution {

    public Input<RealParameter> beta_input = new Input<RealParameter>("beta", "the transmission rate parameter.", Input.Validate.REQUIRED);
    public Input<RealParameter> delta_input = new Input<RealParameter>("delta", "the recovery rate parameter.", Input.Validate.REQUIRED);
    public Input<RealParameter> S0_input = new Input<RealParameter>("S0", "the S_0 at the time of the origin.", Input.Validate.REQUIRED);
    public Input<RealParameter> T_input = new Input<RealParameter>("T", "the origin.", Input.Validate.REQUIRED);

    TreeIntervals intervals;
    VolzSIR volzSIR;

    int integrationStepCount = 1000;

    @Override
    public void initAndValidate() throws Exception {
        intervals = treeIntervals.get();
        if (intervals == null) {
            throw new Exception("Expected treeIntervals to be specified");
        }

        volzSIR = new VolzSIR();
        volzSIR.initByName(
                "n_S0", S0_input.get(),
                "beta", beta_input.get(),
                "gamma", delta_input.get(),
                "origin", T_input.get(),
                "integrationStepCount", integrationStepCount);

        calculateLogP();
    }

    public double calculateLogP() throws Exception {

        FirstOrderIntegrator integrator = new ClassicalRungeKuttaIntegrator(0.001);

        int n = treeIntervals.get().getIntervalCount() + 1;

        A a = new A();

        double[] B = new double[1];


        double logL = 0;

        double[] A = new double[]{n};
        double s = 0.0;
        for (int i = 0; i < intervals.getIntervalCount(); i++) {
            double t = intervals.getInterval(i) + s;

            if (intervals.getIntervalType(i) == IntervalType.COALESCENT) {

                integrator.integrate(a, s, A, t, B);

                // compute derivatives at time t (end of this coalescent interval)
                double[] adot = new double[1];
                a.computeDerivatives(t, B, adot);

                logL += Math.log(-adot[0]);

                A[0] = B[0];
                s = t;
            }
        }

        A[0] = n;

        double T = T_input.get().getValue();

        integrator.integrate(a, 0, A, T, B);

        logL -= (n - 1) * Math.log(n - B[0]);

        return logL;
    }


    class A implements FirstOrderDifferentialEquations {

        @Override
        public int getDimension() {
            return 1; //A (the expected lineage count through time)
        }

        @Override
        public void computeDerivatives(double t, double[] A, double[] Adot)
                throws MaxCountExceededException, DimensionMismatchException {

            double beta = beta_input.get().getValue();

            double S = volzSIR.getNS(t);
            double I = volzSIR.getNI(t);

            // A
            // f_SI = beta * S * I
            Adot[0] = -beta * S * A[0] * A[0] / I;
        }
    }
}
