package phylodynamics.epidemiology;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.TreeDistribution;
import beast.base.evolution.tree.IntervalType;
import beast.base.evolution.tree.TreeIntervals;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.ClassicalRungeKuttaIntegrator;


/**
 * @author Alexei Drummond
 * @author Tim Vaughan
 * @author Alex Popinga
 */

@Description("Calculates the Volz-2009 probability of a beast.tree.")
public class Volz2009TreeDistribution extends TreeDistribution {

    public Input<DeterministicCoalescentSIR> volzSIRInput = new Input<DeterministicCoalescentSIR>("volzSIR",
            "Volz parameters for deterministic coalescent SIR", Input.Validate.REQUIRED);

    TreeIntervals intervals;

    @Override
    public void initAndValidate() {
        intervals = treeIntervalsInput.get();
        if (intervals == null) {
            throw new RuntimeException("Expected treeIntervals to be specified");
        }

        calculateLogP();
    }


    public double calculateLogP() {

        boolean reject = volzSIRInput.get().update();

        // if origin is being sampled it must be constrained to be older than the root of the tree.
        if (volzSIRInput.get().originParameter.get().getValue() < treeIntervalsInput.get().getTotalDuration()) {
            logP = Double.NEGATIVE_INFINITY;
            return logP;
        }

        // if maximal count exceeded or number too small
        if (reject) {
            logP = Double.NEGATIVE_INFINITY;
            return logP;
        }

        FirstOrderIntegrator integrator = new ClassicalRungeKuttaIntegrator(0.01);

        int n = treeIntervalsInput.get().getIntervalCount() + 1;

        A a = new A();

        double[] B = new double[1];
        double logL = 0;
        double[] A = new double[]{n};
        double s = 0.0;


        // Get sampling times and piecewise integration of a(t)
        // Evaluate C=\sum_k(C_k), the coalescent events within interval k

        double C = 0.0;
        for (int i = 0; i < intervals.getIntervalCount(); i++) {
            double t = intervals.getInterval(i) + s;

                if (intervals.getInterval(i) > 0.0) {
                    integrator.integrate(a, s, A, t, B);
                } else {
                    B[0] = A[0];
                }

            C += A[0] - B[0];

            A[0] = B[0];
            s = t;

            if (intervals.getIntervalType(i) == IntervalType.SAMPLE) {
                A[0] += 1;
            }

            // Get coalescent times, evaluate -adot(t) at each coalescent time
            if (intervals.getIntervalType(i) == IntervalType.COALESCENT) {

                double[] adot = new double[1];
                a.computeDerivatives(t, B, adot);

                logL += Math.log(-adot[0]);

            }

        }

        // Final calculations, between tMRCA (last interval) and origin (T)

        double T = volzSIRInput.get().originParameter.get().getValue();
        double u = intervals.getTotalDuration();
        double x = T - u;

        integrator.integrate(a, u, A, T, B);
        C += A[0] - B[0];

        logL -= (n - 1 + x) * Math.log(C);

        logP = logL;
        return logP;


    }

    class A implements FirstOrderDifferentialEquations {

        @Override
        public int getDimension() {
            return 1; //A (the expected lineage count through time)
        }

        @Override
        public void computeDerivatives(double t, double[] A, double[] Adot)
                throws MaxCountExceededException, DimensionMismatchException {

            double beta = volzSIRInput.get().beta();

            double S = volzSIRInput.get().getNS(t);
            double I = volzSIRInput.get().getNI(t);

            // A
            // f_SI = beta * S * I
            Adot[0] = -beta * S * A[0] * A[0] / I;
        }
    }
}
