package phylodynamics.evolution.operators;

import beast.base.core.Description;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.core.Input;
import beast.base.util.Randomizer;
import beast.base.evolution.operator.ScaleOperator;
import beast.base.evolution.tree.Tree;

/**
 * User: Denise
 * Date: Jul 2, 2013
 * Time: 4:03:43 PM
 */

@Description("Scale a parameter that must always be larger than the height of the tree")
public class ScaleOffsetOperator extends ScaleOperator {


    public final Input<Tree> offsetTree = new Input<Tree>("offsetTree", "tree to specify the tree height that is used as offset", Input.Validate.REQUIRED);

    public Input<RealParameter> scalingOffset = new Input<RealParameter>("scalingOffset", "parameter determining the offset z for scaling: x' = z + b * x", Input.Validate.XOR, offsetTree);


    // the offset z for scaling parameter x: x' = z + b * x
    double offset;

    @Override
    public void initAndValidate() {

        super.initAndValidate();

        if (isTreeScaler()) throw new UnsupportedOperationException();

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

        if (offsetTree.get() !=null) offset = offsetTree.get().getRoot().getHeight();

        else offset = scalingOffset.get().getValue();

        try {

            double hastingsRatio;
            final double scale = getScaler();

            // scale a parameter
            final boolean bScaleAll = scaleAllInput.get();
            final boolean bScaleAllIndependently = scaleAllIndependentlyInput.get();

            final RealParameter param = parameterInput.get();

            assert param.getLower() != null && param.getUpper() != null;

            final int dim = param.getDimension();

            if (bScaleAllIndependently) {

                throw new UnsupportedOperationException();

            } else if (bScaleAll) {

                throw new UnsupportedOperationException();

                }
//            } else {
                hastingsRatio = -Math.log(scale);

                // which position to scale
                final int index;
                final BooleanParameter indicators = indicatorInput.get();
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
//            }

            return hastingsRatio;

        } catch (Exception e) {
            // whatever went wrong, we want to abort this operation...
            return Double.NEGATIVE_INFINITY;
        }
    }


}
