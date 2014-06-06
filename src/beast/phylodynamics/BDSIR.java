package beast.phylodynamics;

import beast.core.Citation;
import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.speciation.BirthDeathSkylineModel;
import beast.evolution.tree.TreeInterface;

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
public class BDSIR extends BirthDeathSkylineModel {


    public Input<Function> S0_input =
            new Input<Function>("S0", "The numbers of susceptible individuals");

    public Input<RealParameter> m_dS =
            new Input<RealParameter>("dS", "dS vector containing the changes in numbers of susceptibles per location", Input.Validate.REQUIRED);
    public Input<RealParameter> m_dE =
            new Input<RealParameter>("dE", "dE vector containing the changes in numbers of exposed per location");
    public Input<RealParameter> m_dR =
            new Input<RealParameter>("dR", "dR vector containing the changes in numbers of recovered per location", Input.Validate.REQUIRED);

    public Input<Boolean> checkTreeConsistent = new Input<Boolean>("checkTreeConsistent", "check if trajectory is consistent with number of lineages in tree? default true", true);

    public Input<Boolean> isSeasonal = new Input<Boolean>("isSeasonal", "Is this a SeasonalSIRSEpidemic? default false", false);

    Double S0;
    Double[] dS;
    Double[] dE;
    Double[] dR;

    int dim;
    double T;
    int ntaxa;


    @Override
    public void initAndValidate() throws Exception {

        S0 = (S0_input.get().getArrayValue());

        dS = m_dS.get().getValues();
        dim = dS.length;

        birthChanges = dim - 1;
        super.initAndValidate();

        if (transform) {
            if (R0.get().getDimension() != 1 && !isSeasonal.get())// || becomeUninfectiousRate.get().getDimension() != 1 || samplingProportion.get().getDimension() != 1)
                throw new RuntimeException("R0, becomeUninfectiousRate and samplingProportion have to be 1-dimensional!");
        } else {
            if (birthRate.get().getDimension() != 1 && !isSeasonal.get())//  || death.length != 1 || psi.length != 1)
                throw new RuntimeException("Birth, death and sampling rate have to be 1-dimensional!");
        }

        // todo: add check that intervaltimes make sense (removed for BDSIR in bdsky to allow seasonality)

        T = origin.get().getValue();
        ntaxa = treeInput.get().getLeafNodeCount();

    }


    @Override
    public Double updateRatesAndTimes(TreeInterface tree) {

        super.updateRatesAndTimes(tree);

        T = origin.get().getValue();
        ntaxa = tree.getLeafNodeCount();

        S0 = (S0_input.get().getArrayValue());

        dS = m_dS.get().getValues();

        dE = (m_dE.get() != null) ? m_dE.get().getValues() : (new Double[dS.length]);
        if (dE[0] == null) Arrays.fill(dE, 0.);

        dR = m_dR.get().getValues();


        double cumS = S0 - 1;


        double time;

        Double[] birthSIR = new Double[dim];
        double I = 1.;
        double R = 0.;


        int season = (!isSeasonal.get()) ? 0 : getSeason(T);
        int initialSeason = season;

        if (isSeasonal.get())
            birth[1] = transform ? (R0.get().getValue(1) * becomeUninfectiousRate.get().getValue()) : birthRate.get().getValue(1);

        birthSIR[0] = birth[season] / S0 * cumS;
        for (int i = 0; i < dim - 1; i++) {

            time = (i + 1.) / dim * T;
            if (isSeasonal.get()) season = (initialSeason + getSeason(T - time)) % 2;

            cumS -= dS[i];
            birthSIR[i + 1] = birth[season] / S0 * cumS;

            I += dS[i] - dE[i] - dR[i];
            R += dR[i];

            if (checkTreeConsistent.get() && (I <= 0. || I < lineageCountAtTime(T - time, tree)))
                return Double.NEGATIVE_INFINITY;

        }

        if (cumS < 0 || S0 - cumS < treeInput.get().getLeafNodeCount() || S0 != (cumS + I + R))
            return Double.NEGATIVE_INFINITY;

        adjustBirthRates(birthSIR);
        return 0.;

    }

    /**
     * @param birthSIR
     */
    public void adjustBirthRates(Double[] birthSIR) {

        for (int i = 0; i < totalIntervals; i++) {
            birth[i] = birthSIR[birthChanges > 0 ? index(times[i], birthRateChangeTimes) : 0];
        }
    }


    int getSeason(double time) {   // this assumes that the second minus first change time entry in the xml defines the length of a season

        double seasonLength = birthRateChangeTimesInput.get().getValue(1) - birthRateChangeTimesInput.get().getValue(0);

        double t = (time - birthRateChangeTimesInput.get().getValue(0));

        return (int) Math.floor(1 + t / seasonLength) % 2;

    }


    @Override
    public Boolean isBDSIR() {
        return true;
    }

    public Boolean isSeasonalBDSIR() {
        return isSeasonal.get();
    }


    public int getSIRdimension() {
        return dim;
    }
}
