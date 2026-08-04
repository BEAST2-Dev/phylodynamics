/*
 * Copyright (C) 2019-2025 ETH Zurich
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package phylodynamics.parameterization;

import bdmmprime.parameterization.EpiParameterization;
import bdmmprime.parameterization.SkylineVectorParameter;
import beast.base.inference.parameter.RealParameter;
import phylodynamics.BDSIR;

import java.util.Arrays;

/**
 * Parameterization for BDSIR model, reusing Bdmmprime's EpiParameterization
 */
public class BDSIRParameterization {

    //dirty flag so the values can be refreshed
    public static class DirtySkylineParCheck extends SkylineVectorParameter {
        public DirtySkylineParCheck(RealParameter changeTimes, RealParameter values) {
            super(changeTimes, values);
        }
        public void markDirty() {
			isDirty = true; 
		}
    }

    public static class DirtyEpiParameterization extends EpiParameterization {
        public void markDirty() {
			dirty = true; 
		}
    }

    private final BDSIR bdsir;

    private RealParameter reValues, changeTimes, becomeUninfectious, samplingProp, removalProb;
    private DirtySkylineParCheck reSkyline, buSkyline, spSkyline, rpSkyline;
    private DirtyEpiParameterization epiParameterization;

    private boolean initialised = false;

    public BDSIRParameterization(BDSIR bdsir) {
        this.bdsir = bdsir;
        build(1);
    }

    public DirtyEpiParameterization get() {
        return epiParameterization; //The EpiParameterization for BirthDeathMigrationDistribution.

    }

    private void build(int dim) {
        Double[] re = new Double[dim];
        Arrays.fill(re, 1.0);
        reValues = new RealParameter(re);

        if (dim > 1) {
            Double[] times = new Double[dim - 1];
            for (int i = 0; i < dim - 1; i++) times[i] = (i + 1.0);
            changeTimes = new RealParameter(times);
        } else {
            changeTimes = null;
        }

        becomeUninfectious = new RealParameter(new Double[]{1.0});
        samplingProp = new RealParameter(new Double[]{0.1});
        removalProb = new RealParameter(new Double[]{1.0});   // non-SA model

        reSkyline = new DirtySkylineParCheck(changeTimes, reValues);
        buSkyline = new DirtySkylineParCheck(null, becomeUninfectious);
        spSkyline = new DirtySkylineParCheck(null, samplingProp);
        rpSkyline = new DirtySkylineParCheck(null, removalProb);

        epiParameterization = new DirtyEpiParameterization();
        epiParameterization.ReInput.setValue(reSkyline, epiParameterization);
        epiParameterization.becomeUninfectiousRateInput.setValue(buSkyline, epiParameterization);
        epiParameterization.samplingProportionInput.setValue(spSkyline, epiParameterization);
        epiParameterization.removalProbInput.setValue(rpSkyline, epiParameterization);
    }

    //Called from BDSIR's initAndValidate(), needs dim and T
    public void initialise(int dim, double T) {
        build(dim);

        // changetimes like in bdsky
        if (changeTimes != null) {
            for (int i = 0; i < dim - 1; i++)
                changeTimes.setValue(i, (i + 1.) / dim * T);
        }

        epiParameterization.processLengthInput.setValue(bdsir.origin.get(), epiParameterization);
        epiParameterization.initAndValidate();
        initialised = true;
    }


    public void SIRValues() {

        Double result = bdsir.updateRatesAndTimes(bdsir.treeInput.get());
        bdsir.treeConsistent = (result != Double.NEGATIVE_INFINITY);

        double b = bdsir.becomeUninfectiousRate.get().getArrayValue();
        double p = bdsir.samplingProportion.get().getArrayValue();

        for (int i = 0; i < bdsir.dim; i++)
            reValues.setValue(i, bdsir.birth[i] / b);

        if (changeTimes != null) {
            for (int i = 0; i < bdsir.dim - 1; i++)
                changeTimes.setValue(i, (i + 1.) / bdsir.dim * bdsir.T);
        }

        becomeUninfectious.setValue(0, b);
        samplingProp.setValue(0, p);

        reSkyline.markDirty();
        buSkyline.markDirty();
        spSkyline.markDirty();
        rpSkyline.markDirty();
        epiParameterization.markDirty();
    }
}
