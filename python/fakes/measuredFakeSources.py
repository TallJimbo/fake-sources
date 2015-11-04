
import lsst.afw.image
import lsst.afw.detection
import lsst.pex.config
from lsst.pipe.tasks.fakes import FakeSourcesConfig, FakeSourcesTask


class MeasuredFakeSourcesConfig(FakeSourcesConfig):
    measurement = lsst.pex.config.ConfigurableField(
        target=SourceMeasurementTask,
        doc=("Subtask to run measurement algorithms and blendedness metrics on "
             "the isolated fake sources before they're added to the full image.")
    )

    def setDefaults(self):
        FakeSourcesConfig.setDefaults(self)
        self.measurement.blendedness.doOld = False
        self.measurement.doBlendedness = True
        self.measurement.algorithms.names = ["shape.hsm.moments"]  # could use shape.sdss if HSM unavailable
        self.measurement.slots.psfFlux = None
        self.measurement.slots.apFlux = None
        self.measurement.slots.modelFlux = None
        self.measurement.slots.instFlux = None
        self.measurement.slots.shape = "shape.hsm.moments"


class MeasuredFakeSourcesTask(FakeSourcesTask):
    """An intermediate base class for fake sources tasks that provides methods to perform "ideal"
    measurements on fake sources as they're added to images, including blendedness metrics.

    A subclass of MeasuredFakeSourcesTask should implement its run() method as usual, but call
    methods of its base class at three points:
     - call makeCatalog before starting a loop over all sources to be added to an image;
     - call measureFake as each object is added, with a postage stamp containing just the fake object;
     - call measureParents after all objects have been added to the image.
    """

    ConfigClass = MeasuredFakeSourcesConfig

    def __init__(self, **kwargs):
        FakeSourcesTask.__init__(self, **kwargs)
        self.schema = lsst.afw.table.SourceTable.makeMinimalSchema()
        self.algMetadata = PropertyList()
        self.makeSubTask("measurement", schema=self.schema, algMetadata=self.algMetadata)


    def makeCatalog(self):
        """Create a catalog to hold the results of measuring fake objects in ideal postage stamps.

        Should be called by fake sources tasks before starting a loop over objects.  The returned
        catalog should then be passed to measureFake() as each object is added, and finally to
        measureParents() after all sources are added.
        """
        catalog = lsst.afw.table.SourceCatalog(self.schema)
        self.measurement.config.slots.setupTable(catalog.table, prefix=self.config.prefix)
        return catalog

    def measureFake(self, realExposure, fakeImage, xyPosition, catalog):
        """Run measurement algorithms on a single fake object in its isolated postage stamp.

        This includes running the portion of the blendedness metrics that operate on the child
        image.

        Because this runs on an ideal image with no noise (and accordingly zeros in the variance plane)
        some algorithms may not perform well, particularly in the calculation of their uncertainties.

        Returns a SourceRecord containing the measurements.
        """
        # Make a fake Exposure out of the postage stamp with the metadata from the real one.
        fakeMaskedImage = lsst.afw.image.MaskedImageF(fakeImage)
        fakeExposure = realExposure.Factory(fakeMaskedImage)
        fakeExposure.setPsf(realExposure.getPsf())
        fakeExposure.setWcs(realExposure.getWcs())
        fakeExposure.setCalib(realExposure.getCalib())
        fakeExposure.setFilter(realExposure.getFilter())
        # Create a new record to save the idealized measurements in.
        record = catalog.addNew()
        # We just use the bounding box of the postage stamp as the Footprint for now; could
        # be more sophisticated we were running algorithms that required that.
        footprint = lsst.afw.detection.Footprint(fakeImage.getBBox(lsst.afw.image.PARENT))
        record.setFootprint(footprint)
        # Run measurement algorithms.  This is a bit fragile, as we're depending a bit on some
        # implementation details of SourceMeasurementTask, but it should work fine for now.
        self.measurement.measurer.apply(record, fakeExposure, xyPosition)
        # Run blendedness on the child; measuring the parent will have to wait for later.
        if self.config.doBlendedness:
            self.measurement.blendedness.measureChildPixels(fakeMaskedImage, record)
        return record

    def measureParents(self, realExposure, catalog):
        """After adding all fake sources to an image, finish blendedness measurements.

        Should be called by fake source tasks after adding all sources to the real image, and calling
        measureFake on each of them as they're added.
        """
        if self.config.doBlendedness
            for record in catalog:
                self.measurement.blendedness.measureParentPixels(realExposure.getMaskedImage(), record)
