package beast.evolution.operators;

import sun.reflect.generics.reflectiveObjects.NotImplementedException;
import beast.core.parameter.RealParameter;
import beast.core.parameter.BooleanParameter;
import beast.core.Input;
import beast.util.Randomizer;
import beast.evolution.tree.Tree;

/**
 * User: Denise
 * Date: Jul 2, 2013
 * Time: 4:03:43 PM
 */
public class ScaleOffsetOperator extends ScaleOperator{


    public final Input<Tree> offsetTree = new Input<Tree>("offsetTree", "tree to specify the tree height that is used as offset", Input.Validate.REQUIRED);

    public Input<RealParameter> scalingOffset = new Input<RealParameter>("scalingOffset", "parameter determining the offset z for scaling: x' = z + b * x", Input.Validate.XOR, offsetTree);


    // the offset z for scaling parameter x: x' = z + b * x
    double offset;

    @Override
    public void initAndValidate() throws Exception {

        super.initAndValidate();

        if (m_bIsTreeScaler) throw new NotImplementedException();

        if (offsetTree.get() !=null) offset = offsetTree.get().getRoot().getHeight();

        else offset = scalingOffset.get().getValue();
    }


    /**
     * override this for proposals,
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {
        try {

            double hastingsRatio;
            final double scale = getScaler();

            // scale a parameter
            final boolean bScaleAll = m_pScaleAll.get();
            final int nDegreesOfFreedom = m_pDegreesOfFreedom.get();
            final boolean bScaleAllIndependently = m_pScaleAllIndependently.get();

            final RealParameter param = m_pParameter.get(this);

            assert param.getLower() != null && param.getUpper() != null;

            final int dim = param.getDimension();

            if (bScaleAllIndependently) {
                // update all dimensions independently.
                hastingsRatio = 0;
                for (int i = 0; i < dim; i++) {

                    final double scaleOne = getScaler();
                    final double newValue = offset + scaleOne * (param.getValue(i) - offset);

                    hastingsRatio -= Math.log(scaleOne);

                    if (outsideBounds(newValue, param)) {
                        return Double.NEGATIVE_INFINITY;
                    }

                    param.setValue(i, newValue);
                }
            } else if (bScaleAll) {
                // update all dimensions
                // hasting ratio is dim-2 times of 1dim case. would be nice to have a reference here
                // for the proof. It is supposed to be somewhere in an Alexei/Nicholes article.
                final int df = (nDegreesOfFreedom > 0) ? -nDegreesOfFreedom : dim - 2;
                hastingsRatio = df * Math.log(scale);

                // all Values assumed independent!
                for (int i = 0; i < dim; i++) {
                    final double newValue = offset + scale * (param.getValue(i) - offset);

                    if (outsideBounds(newValue, param)) {
                        return Double.NEGATIVE_INFINITY;
                    }
                    param.setValue(i, newValue);
                }
            } else {
                hastingsRatio = -Math.log(scale);

                // which position to scale
                final int index;
                final BooleanParameter indicators = m_indicator.get();
                if (indicators != null) {
                    final int nDim = indicators.getDimension();
                    Boolean[] indicator = indicators.getValues();
                    final boolean impliedOne = nDim == (dim - 1);

                    // available bit locations. there can be hundreds of them. scan list only once.
                    final int[] loc = new int[nDim + 1];
                    int nLoc = 0;

                    if (impliedOne) {
                        loc[nLoc] = 0;
                        ++nLoc;
                    }
                    for (int i = 0; i < nDim; i++) {
                        if (indicator[i]) {
                            loc[nLoc] = i + (impliedOne ? 1 : 0);
                            ++nLoc;
                        }
                    }

                    if (nLoc > 0) {
                        final int rand = Randomizer.nextInt(nLoc);
                        index = loc[rand];
                    } else {
                        return Double.NEGATIVE_INFINITY; // no active indicators
                    }

                } else {
                    // any is good
                    index = Randomizer.nextInt(dim);
                }

                final double oldValue = param.getValue(index);

                if (oldValue == 0) {
                    // Error: parameter has value 0 and cannot be scaled
                    return Double.NEGATIVE_INFINITY;
                }

                final double newValue = offset + scale * (oldValue - offset);

                if (outsideBounds(newValue, param)) {
                    // reject out of bounds scales
                    return Double.NEGATIVE_INFINITY;
                }

                param.setValue(index, newValue);
                // provides a hook for subclasses
                //cleanupOperation(newValue, oldValue);
            }

            return hastingsRatio;

        } catch (Exception e) {
            // whatever went wrong, we want to abort this operation...
            return Double.NEGATIVE_INFINITY;
        }
    }


}
