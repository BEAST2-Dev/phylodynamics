package phylodynamics;

import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import bdmmprime.distribution.BirthDeathMigrationDistribution;
import phylodynamics.parameterization.BDSIRParameterization;
import beast.base.evolution.tree.TreeInterface;
import org.apache.commons.math.special.Gamma;

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

    public Input<Function> origin =
            new Input<Function>("origin", "The origin of infection", Input.Validate.REQUIRED);	// Previously inherited from bdsky's BirthDeathSkylineModel
    public Input<Function> reproductiveNumberInput =
            new Input<Function>("reproductiveNumber", "The basic reproduction number", Input.Validate.REQUIRED); 	// Previously inherited from bdsky's BirthDeathSkylineModel
    public Input<Function> becomeUninfectiousRate =
            new Input<Function>("becomeUninfectiousRate", "Rate at which individuals become uninfectious (death + sampling)", Input.Validate.REQUIRED); 	// Previously inherited from bdsky's BirthDeathSkylineModel
    public Input<Function> samplingProportion =
            new Input<Function>("samplingProportion", "Proportion of samples taken at the time of become uninfectious", Input.Validate.REQUIRED);		// Previously inherited from bdsky's BirthDeathSkylineModel

    Double S0;
    Double[] dS;
    Double[] dE;
    Double[] dR;

    public int dim;
    public double T;
    int ntaxa;

    public double[] birth;
    int birthChanges;
    public boolean treeConsistent = true;

	// BDSIR's internal parameterization that supplies birth/death/sampling rates to BirthDeathMigrationDistribution.
    private final BDSIRParameterization bdsirParameterization = new BDSIRParameterization(this);

    public BDSIR() {
        super();
        parameterizationInput.setValue(bdsirParameterization.get(), this);
    }

    @Override
    public void initAndValidate() {

        S0 = (S0_input.get().getArrayValue());

        dS = m_dS.get().getValues();
        dim = dS.length;

        birthChanges = dim - 1;
		birth = new double[dim];

        T = origin.get().getArrayValue();
        bdsirParameterization.initialise(dim, T);
        parameterizationInput.setValue(bdsirParameterization.get(), this);
        super.initAndValidate();

        /*if (transform) {
            if (reproductiveNumberInput.get().getDimension() != 1 && !isSeasonal.get())// || becomeUninfectiousRate.get().getDimension() != 1 || samplingProportion.get().getDimension() != 1)
                throw new RuntimeException("R0, becomeUninfectiousRate and samplingProportion have to be 1-dimensional!");
        } else {
            if (birthRate.get().getDimension() != 1 && !isSeasonal.get())//  || death.length != 1 || psi.length != 1)
                throw new RuntimeException("Birth, death and sampling rate have to be 1-dimensional!");
        }*/

        // todo: add check that intervaltimes make sense (removed for BDSIR in bdsky to allow seasonality)

        T = origin.get().getArrayValue();
        ntaxa = treeInput.get().getLeafNodeCount();

    }


    //@Override
    public Double updateRatesAndTimes(TreeInterface tree) {

        // super.updateRatesAndTimes(tree); no such method in BDMM Prime

        T = origin.get().getArrayValue();
        ntaxa = tree.getLeafNodeCount();

        S0 = (S0_input.get().getArrayValue());

        dS = m_dS.get().getValues();

        dE = (m_dE.get() != null) ? m_dE.get().getValues() : (new Double[dS.length]);
        if (dE[0] == null) Arrays.fill(dE, 0.);

        dR = m_dR.get().getValues();


        double cumS = S0 - 1;


        double time;

        double[] birthSIR = new double[dim];
        double I = 1.;
        double R = 0.;


        int season = (!isSeasonal.get()) ? 0 : getSeason(T);
        int initialSeason = season;

		// needed for season = 0(default), as it will cause -Infinity, was earlier in bdsky
        birth[0] = reproductiveNumberInput.get().getArrayValue(0) * becomeUninfectiousRate.get().getArrayValue();
        if (isSeasonal.get())
            // birth[1] = transform ? (reproductiveNumberInput.get().getArrayValue(1) * becomeUninfectiousRate.get().getArrayValue()) : birthRate.get().getArrayValue(1);
			birth[1] = reproductiveNumberInput.get().getArrayValue(1) * becomeUninfectiousRate.get().getArrayValue();

        birthSIR[0] = birth[season] / S0 * cumS;
        for (int i = 0; i < dim - 1; i++) {

            time = (i + 1.) / dim * T;
            if (isSeasonal.get()) season = (initialSeason + getSeason(T - time)) % 2;

            cumS -= dS[i];
            birthSIR[i + 1] = birth[season] / S0 * cumS;

            I += dS[i] - dE[i] - dR[i];
            R += dR[i];

			// lineageCountAtTime() was earlier inherited from bdsky. Same method is declared directly here in BDSIR.
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
    public void adjustBirthRates(double[] birthSIR) {

        /*for (int i = 0; i < totalIntervals; i++) {
            birth[i] = birthSIR[birthChanges > 0 ? index(times[i], birthRateChangeTimes) : 0];*/
        for (int i = 0; i < dim; i++) {
            birth[i] = birthSIR[i];
        }
    }

    @Override
    public double calculateTreeLogLikelihood(TreeInterface tree) {
		// recompute the SIR trajectory before BDMM Prime runs
        bdsirParameterization.SIRValues();
        if (!treeConsistent)
            return Double.NEGATIVE_INFINITY;

        double logP = super.calculateTreeLogLikelihood(tree);
        int internalNodeCount = tree.getLeafNodeCount() - ((beast.base.evolution.tree.Tree) tree).getDirectAncestorNodeCount() - 1;
        logP -= Math.log(2) * internalNodeCount;
        logP += Gamma.logGamma(tree.getLeafNodeCount() + 1);
        return logP;
    }

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

        //double seasonLength = birthRateChangeTimesInput.get().getValue(1) - birthRateChangeTimesInput.get().getValue(0);

        //double t = (time - birthRateChangeTimesInput.get().getValue(0));

        double seasonLength = T / dim;
        double t = time;
        return (int) Math.floor(1 + t / seasonLength) % 2;

    }


    //@Override
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
