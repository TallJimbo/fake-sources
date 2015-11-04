import os
import fcntl
import numpy as np

import lsst.afw.image
import lsst.afw.geom
import lsst.afw.math
import lsst.afw.cameraGeom
import lsst.afw.coord
import lsst.pex.config
import pyfits as fits

import makeFakeGalaxy as makeFake
import FakeSourceLib as fsl

from .measuredFakeSources import MeasuredFakeSourcesTask, MeasuredFakeSourcesConfig

class PositionGalSimFakesConfig(MeasuredFakeSourcesConfig):
    galList = lsst.pex.config.Field(dtype=str, doc="catalog of galaxies to add")
    maxMargin = lsst.pex.config.Field(dtype=int, default=600, optional=True,
                                   doc="Size of margin = 1/2 of maximum galsim image size, in pixels")
    seed = lsst.pex.config.Field(dtype=int, default=1,
                                 doc="Seed for random number generator")
    galType = lsst.pex.config.ChoiceField(dtype=str, default='sersic',
                                          allowed={'dsersic':'double sersic galaxies added',
                                                   'sersic':'single sersic galaxies added',
                                                   'real':'real HST galaxy images added'},
                                          doc='type of GalSim galaxies to add')
    addShear = lsst.pex.config.Field(dtype=bool, default=False,
                                     doc='include shear in the galaxies')


class PositionGalSimFakesTask(MeasuredFakeSourcesTask):
    ConfigClass = PositionGalSimFakesConfig

    def __init__(self, **kwargs):
        MeasuredFakeSourcesTask.__init__(self, **kwargs)
        print "RNG seed:", self.config.seed
        self.rng = lsst.afw.math.Random(self.config.seed)
        self.npRand = np.random.RandomState(self.config.seed)
        self.galData = fits.open(self.config.galList)[1].data


    def run(self, exposure, background):

        self.log.info("Adding fake galaxies at real positions")
        psf = exposure.getPsf()
        md = exposure.getMetadata()
        expBBox = exposure.getBBox(lsst.afw.image.PARENT)
        wcs = exposure.getWcs()
        skyToPixelMatrix = wcs.getLinearTransform().invert().getMatrix() / 3600.0

        """ Deal with the skipped ones """
        skipLog = 'runAddFake.skipped'
        if not os.path.isfile(skipLog):
            dum = os.system('touch ' + skipLog)

        measCat = self.makeCatalog()

        for igal, gal in enumerate(self.galData):
            try:
                galident = gal["ID"]
            except KeyError:
                galident = igal + 1

            try:
                flux = exposure.getCalib().getFlux(float(gal['mag']))
            except KeyError:
                raise KeyError("No mag column in %s table"%self.config.galList)

            try:
                galCoord = lsst.afw.coord.Coord(lsst.afw.geom.Point2D(gal['RA'], gal['DEC']))
            except KeyError:
                raise("No RA/DEC column in %s table"%self.config.galList)

            galXY = wcs.skyToPixel(galCoord)
            bboxI = exposure.getBBox(lsst.afw.image.PARENT)
            bboxI.grow(self.config.maxMargin)
            if not bboxI.contains(lsst.afw.geom.Point2I(galXY)):
                #self.log.info("Not Useful: Skipping fake %d" % galident)
                continue

            #this is extrapolating for the PSF, probably not a good idea
            psfImage = psf.computeKernelImage(galXY)
            try:
                galArray = makeFake.makeGalaxy(flux, gal, psfImage.getArray(),
                                               self.config.galType,
                                               transform = skyToPixelMatrix,
                                               addShear=self.config.addShear)
            except KeyError as kerr:
                self.log.info("GalSim Key Error: Skipping fake %d"%galident)
                self.log.info("       Mag, nSer, Reff, b/a: %5.2f, %5.2f, %6.2f, %4.1f"%(
                    gal['mag'], gal['sersic_n'], gal['reff'], gal['b_a']))
                self.log.info(kerr.message)
                with open(skipLog, "a") as slog:
                    try:
                        fcntl.flock(slog, fcntl.LOCK_EX)
                        slog.write("%8d , galsimK\n"%galident)
                        fcntl.flock(slog, fcntl.LOCK_UN)
                    except IOError:
                        continue
                continue
            except ValueError as verr:
                self.log.info("GalSim Value Error: Skipping fake %d"%galident)
                self.log.info("       Mag, nSer, Reff, b/a: %5.2f, %5.2f, %6.2f, %4.1f"%(
                    gal['mag'], gal['sersic_n'], gal['reff'], gal['b_a']))
                self.log.info(verr.message)
                with open(skipLog, "a") as slog:
                    try:
                        fcntl.flock(slog, fcntl.LOCK_EX)
                        slog.write("%8d , galsimV\n"%galident)
                        fcntl.flock(slog, fcntl.LOCK_UN)
                    except IOError:
                        continue
                continue
            except RuntimeError as rerr:
                self.log.info("GalSim Runtime Error: Skipping fake %d"%galident)
                self.log.info("       Mag, nSer, Reff, b/a: %5.2f, %5.2f, %6.2f, %4.1f"%(
                    gal['mag'], gal['sersic_n'], gal['reff'], gal['b_a']))
                self.log.info(rerr.message)
                with open(skipLog, "a") as slog:
                    try:
                        fcntl.flock(slog, fcntl.LOCK_EX)
                        slog.write("%8d , galsimR\n"%galident)
                        fcntl.flock(slog, fcntl.LOCK_UN)
                    except IOError:
                        continue
                continue

            galImage = lsst.afw.image.ImageF(galArray.astype(np.float32))
            galBBox = galImage.getBBox(lsst.afw.image.PARENT)
            galImage = lsst.afw.math.offsetImage(galImage,
                                                 galXY.getX() - galBBox.getWidth()/2.0 + 0.5,
                                                 galXY.getY() - galBBox.getHeight()/2.0 + 0.5,
                                                 'lanczos3')
            galBBox = galImage.getBBox(lsst.afw.image.PARENT)


            #check that we're within the larger exposure, otherwise crop
            if expBBox.contains(galImage.getBBox(lsst.afw.image.PARENT)) is False:
                newBBox = galImage.getBBox(lsst.afw.image.PARENT)
                newBBox.clip(expBBox)
                if newBBox.getArea() <= 0:
                    self.log.info("BBoxEdge Error: Skipping fake %d"%galident)
                    with open(skipLog, "a") as slog:
                        try:
                            fcntl.flock(slog, fcntl.LOCK_EX)
                            slog.write("%8d , bboxEdge\n"%galident)
                            fcntl.flock(slog, fcntl.LOCK_UN)
                        except IOError:
                            continue
                    continue
                self.log.info("Cropping FAKE%d from %s to %s"%(galident, str(galBBox), str(newBBox)))
                galImage = galImage.Factory(galImage, newBBox, lsst.afw.image.PARENT)
                galBBox = newBBox

            self.measureFake(exposure, galImage, galXY, measCat)

            galMaskedImage = fsl.addNoise(galImage, exposure.getDetector(),
                                          rand_gen=self.npRand)

            # Put information of the added fake galaxies into the header
            md.set("FAKE%s" % str(galident), "%.3f, %.3f" % (galXY.getX(), galXY.getY()))
            self.log.info("Adding fake %s at: %.1f,%.1f"% (str(galident), galXY.getX(), galXY.getY()))

            galMaskedImage.getMask().set(self.bitmask)
            subMaskedImage = exposure.getMaskedImage().Factory(exposure.getMaskedImage(),
                                                               galMaskedImage.getBBox(lsst.afw.image.PARENT),
                                                               lsst.afw.image.PARENT)
            subMaskedImage += galMaskedImage

        self.measureParents(exposure, measCat)

        return measCat
