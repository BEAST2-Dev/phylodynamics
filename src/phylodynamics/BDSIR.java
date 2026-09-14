package phylodynamics;

import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import bdmmprime.distribution.BirthDeathMigrationDistribution;
import phylodynamics.parameterization.BDSIRParameterization;
import beast.base.evolution.tree.TreeInterface;

import java.util.Arrays;

/**
 * @author Denise Kuhnert
 */

@Description("Phylodynamic tree prior that couples compartmental models (e.g. SIR, SIRS) " +
        "with a piecewise constant birth-death-sampling process.")
@Citation("Simultaneous reconstruction of evolutionary history and epidemiological "
        + "dynamics from viral sequences with the birth–death SIR model. "
        + "Denise Kuehnert, Tanja Stadler, Timothy Vaughan, and Alexei Drummond, "
        + "J. R. Soc. Interface, 11:20131106 (2014). ")
public class BDSIR extends BirthDeathMigrationDistribution {

    // bdsky's equivalent origin, reproductiveNumber, becomeUnifectiousRate and samplingProportion now declared in BDSIRParameterization.

    Double S0;
    Double[] dS;
    Double[] dE;
    Double[] dR;

    public int dim;
    public double T;
    int ntaxa;

    public double[] birth;
    public double[] times;
    public int totalIntervals;
    public double[] birthSIR;
    int birthChanges;
    public boolean treeConsistent = true;

    private BDSIRParameterization bdsirParameterization;

    @Override
    public void initAndValidate() {

        bdsirParameterization = (BDSIRParameterization) parameterizationInput.get(); // To make S0_input, etc. reachable. Parameterization is otherwise private in BirthDeathMigrationDistribution.
        bdsirParameterization.setBDSIR(this); // Parameterization can now see BDSIR

        S0 = (bdsirParameterization.S0_input.get().getArrayValue());

        dS = bdsirParameterization.m_dS.get().getValues();
        dim = dS.length;

        birthChanges = dim - 1;
        birth = new double[dim];

        // T for building interval boundaries.bdsky's origin.
        T = bdsirParameterization.processLengthInput.get().getArrayValue();

        super.initAndValidate();

        /*if (bdsirParameterization.ReInput.get() != null && bdsirParameterization.becomeUninfectiousRateInput.get() != null && bdsirParameterization.samplingProportionInput.get() != null) {
            if (bdsirParameterization.ReInput. != 1 && !bdsirParameterization.isSeasonal.get())// || becomeUninfectiousRate.get().getDimension() != 1 || samplingProportion.get().getDimension() != 1)
            {
                throw new RuntimeException("R0, becomeUninfectiousRate and samplingProportion have to be 1-dimensional!");
            } else {
                if (bdsirParameterization.getBirthRates().length != 1 && !bdsirParameterization.isSeasonal.get())//  || death.length != 1 || psi.length != 1)
                    throw new RuntimeException("Birth, death and sampling rate have to be 1-dimensional!");
            }
        }*/

        // todo: add check that intervaltimes make sense (removed for BDSIR in bdsky to allow seasonality)

        T = bdsirParameterization.processLengthInput.get().getArrayValue(); //
        ntaxa = treeInput.get().getLeafNodeCount();
    }


    public Double updateRatesAndTimes(TreeInterface tree) {

        T = bdsirParameterization.processLengthInput.get().getArrayValue();
        ntaxa = tree.getLeafNodeCount();

        S0 = (bdsirParameterization.S0_input.get().getArrayValue());

        dS = bdsirParameterization.m_dS.get().getValues();
        times = bdsirParameterization.getIntervalEndTimes();
        totalIntervals = bdsirParameterization.getTotalIntervalCount();
        dim = dS.length;
        birthChanges = dim -1;
        dE = (bdsirParameterization.m_dE.get() != null) ? bdsirParameterization.m_dE.get().getValues() : (new Double[dS.length]);
        if (dE[0] == null) Arrays.fill(dE, 0.);

        dR = bdsirParameterization.m_dR.get().getValues();

        double cumS = S0 - 1;

        double time;

        birthSIR = new double[dim];
        double I = 1.;
        double R = 0.;

        int season = (!bdsirParameterization.isSeasonal.get()) ? 0 : getSeason(T);
        int initialSeason = season;

        birth[0] = bdsirParameterization.ReInput.get().getValuesAtTime(0)[0] * bdsirParameterization.becomeUninfectiousRateInput.get().getValuesAtTime(0)[0];
        if (bdsirParameterization.isSeasonal.get())
            birth[1] = bdsirParameterization.ReInput.get().getValuesAtTime(0)[1] * bdsirParameterization.becomeUninfectiousRateInput.get().getValuesAtTime(0)[1];

        birthSIR[0] = birth[season] / S0 * cumS;
        for (int i = 0; i < dim - 1; i++) {

            time = (i + 1.) / dim * T;
            if (bdsirParameterization.isSeasonal.get()) season = (initialSeason + getSeason(T - time)) % 2;

            cumS -= dS[i];
            birthSIR[i + 1] = birth[season] / S0 * cumS;

            I += dS[i] - dE[i] - dR[i];
            R += dR[i];

            if (bdsirParameterization.checkTreeConsistent.get() && (I <= 0. || I < lineageCountAtTime(T - time, tree)))
                return Double.NEGATIVE_INFINITY;

        }

        if (cumS < 0 || S0 - cumS < treeInput.get().getLeafNodeCount() || S0 != (cumS + I + R))
            return Double.NEGATIVE_INFINITY;

        birth = new double[totalIntervals]; // resize to merged Intervals
        adjustBirthRates(birthSIR);
        return 0.;

    }

    /**
     * @param birthSIR
     */
    public void adjustBirthRates(double[] birthSIR) {
        for (int i = 0; i < totalIntervals; i++) {
            birth[i] = birthSIR[birthChanges > 0 ? index(times[i], bdsirParameterization.ReInput.get().getChangeTimes()) : 0];
        }
    }

    /* To find which interval does t fall into. The method was earlier present in bdsky */
    public int index(double t, double[] times) {
        int epoch = java.util.Arrays.binarySearch(times, t);

        if (epoch < 0)
            epoch = -epoch - 1;
        return epoch;
    }

    @Override
    public double calculateTreeLogLikelihood(TreeInterface tree) {
        return super.calculateTreeLogLikelihood(tree);
    }


    // The method was earlier present in bdsky
    public int lineageCountAtTime(double time, TreeInterface tree) {
        int count = 1;
        int tipCount = tree.getLeafNodeCount();
        for (int i = tipCount; i < tipCount + tree.getInternalNodeCount(); i++) {
            if (tree.getNode(i).getHeight() > time) count += 1;
        }
        for (int i = 0; i < tipCount; i++) {
            if (tree.getNode(i).getHeight() >= time) count -= 1;
        }
        return count;
    }

    int getSeason(double time) {   // this assumes that the second minus first change time entry in the xml defines the length of a season

        double seasonLength = bdsirParameterization.getBirthRateChangeTimes()[1] - bdsirParameterization.getBirthRateChangeTimes()[0];
        double t = (time - bdsirParameterization.getBirthRateChangeTimes()[0]);
        return (int) Math.floor(1 + t / seasonLength) % 2;

    }


    public Boolean isSeasonalBDSIR() {
        return bdsirParameterization.isSeasonal.get();
    }

}
