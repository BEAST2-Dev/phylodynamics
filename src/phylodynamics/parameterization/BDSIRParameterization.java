package phylodynamics.parameterization;

import bdmmprime.parameterization.EpiParameterization;
import phylodynamics.BDSIR;

public class BDSIRParameterization extends EpiParameterization{

    private BDSIR bdsir;

    public void setBDSIR(BDSIR bdsir_) {
        this.bdsir = bdsir_;
    }

    @Override
    public void initAndValidate() {
        if (bdsir == null)
            return;
        super.initAndValidate();
    }


    @Override
    public double[] getBirthRateChangeTimes() {
        Double result = bdsir.updateRatesAndTimes(bdsir.treeInput.get());
        bdsir.treeConsistent = (result != Double.NEGATIVE_INFINITY);

        // same as bdsky‘s getChangeTimes() -> equidistant intervals over the processLength
        int numChanges = Math.max(bdsir.dim-1 , 0);
        double intervalWidth = bdsir.T/bdsir.dim;
        double[] changeTimes = new double[numChanges];
        for (int i =0; i < numChanges; i++)
            changeTimes[i] = intervalWidth*(i+1);
        return changeTimes;
    }

    @Override
    public double[] getBirthRateValues(double time) {
        double intervalWidth = bdsir.T/bdsir.dim;
        int i = 0;
        while (i < bdsir.dim-1 && time > intervalWidth*(i+1))
            i++;
        return new double[]{bdsir.birth[i]};
    }

    @Override
    protected boolean requiresRecalculation(){
        dirty = true;
        return true;
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