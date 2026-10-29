"""Bounded native-resolution RGB reads, preserving the release pixel frame."""
import numpy as np
from PIL import Image


class RGBReader:
    def __init__(self, path, shape_yx):
        self.path = str(path)
        self.height, self.width = map(int, shape_yx)
        self.slide = self.array = self.tf = self.store = None
        try:
            import openslide
            self.slide = openslide.OpenSlide(self.path)
            if self.slide.dimensions != (self.width, self.height):
                self.slide.close()
                self.slide = None
            else:
                self.slide.set_cache(openslide.OpenSlideCache(64 * 1024 ** 2))
                self.backend = 'openslide'
        except Exception as error:
            if error.__class__.__name__ not in ['OpenSlideUnsupportedFormatError', 'OpenSlideError', 'ModuleNotFoundError', 'ImportError']:
                raise
            self.slide = None
        if self.slide is None:
            import tifffile
            self.tf = tifffile.TiffFile(self.path)
            self.axes = self.tf.series[0].axes
            try:
                self.array = tifffile.memmap(self.path, mode='r')
                self.backend = 'tifffile.memmap'
            except (ValueError, TypeError):
                import zarr
                self.store = self.tf.series[0].aszarr()
                self.array = zarr.open(self.store, mode='r')
                self.backend = 'tifffile.zarr'
            assert 'Y' in self.axes and 'X' in self.axes
            assert self.array.shape[self.axes.index('Y')] == self.height
            assert self.array.shape[self.axes.index('X')] == self.width
            assert self.array.dtype == np.uint8, 'H&E reader requires uint8 RGB input'

    def region(self, x0, y0, x1, y1):
        if self.slide is not None:
            rgba = self.slide.read_region((x0, y0), 0, (x1 - x0, y1 - y0))
            background = Image.new('RGB', rgba.size, 'white')
            background.paste(rgba, mask=rgba.getchannel('A'))
            return np.asarray(background)
        key = [0] * len(self.axes)
        key[self.axes.index('Y')] = slice(y0, y1)
        key[self.axes.index('X')] = slice(x0, x1)
        channel = 'S' if 'S' in self.axes else 'C'
        assert channel in self.axes
        assert self.array.shape[self.axes.index(channel)] >= 3
        key[self.axes.index(channel)] = slice(0, 3)
        remaining = ''.join(axis for axis, value in zip(self.axes, key) if isinstance(value, slice))
        data = np.asarray(self.array[tuple(key)])
        data = np.transpose(data, [remaining.index(axis) for axis in ['Y', 'X', channel]])
        assert data.shape == (y1 - y0, x1 - x0, 3)
        return data

    def crop(self, x, y, size=150):
        left, top = int(round(float(x))) - size // 2, int(round(float(y))) - size // 2
        right, bottom = left + size, top + size
        x0, y0 = max(0, left), max(0, top)
        x1, y1 = min(self.width, right), min(self.height, bottom)
        out = np.full((size, size, 3), 255, dtype=np.uint8)
        if x1 > x0 and y1 > y0:
            out[y0 - top:y1 - top, x0 - left:x1 - left] = self.region(x0, y0, x1, y1)
        return Image.fromarray(out, 'RGB')

    def close(self):
        if self.slide is not None:
            self.slide.close()
        if self.store is not None:
            self.store.close()
        if self.tf is not None:
            self.tf.close()
