package phylodynamics.parameterization;

import bdmmprime.parameterization.EpiParameterization;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import phylodynamics.BDSIR;
import java.util.Arrays;

public class BDSIRParameterization extends EpiParameterization{
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

    private BDSIR bdsir; // set using setBDSIR method. getBirthRateValues needs bdsir.birth

    public void setBDSIR(BDSIR bdsir_) {
        this.bdsir = bdsir_;
        requiresRecalculation();
    }

    @Override
    public void initAndValidate() {
        intervalEndTimes = null; // protected in bdmmprime.parameterization
        birthRates = null; // protected in bdmmprime.parameterization
        super.initAndValidate();
    }


    @Override
    public double[] getBirthRateValues(double time) {
        if(bdsir == null)
            return ZERO_VALUE_ARRAY;

        if (time == intervalEndTimes[0]) { // to avoid updateRatesAndTimes running intervalEndTimes.len times
            Double result = bdsir.updateRatesAndTimes(bdsir.treeInput.get());
        }
        return new double[]{bdsir.birth[bdsir.index(time, bdsir.times)]};
    }

    @Override
    public double[] getIntervalEndTimes() {
        return intervalEndTimes;
    }

    @Override
    public int getTotalIntervalCount() {
        return intervalEndTimes.length;
    }


    @Override
    public boolean valuesAreValid(){
        if(super.valuesAreValid() && bdsir.treeConsistent) {
            return true;
        }
        else{
            return false;
        }
    }
}