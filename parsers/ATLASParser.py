import os
import numpy as np
from datetime import datetime, timedelta
import kbmod.search as kb

from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.units as u
from astropy.io import fits

from kbmod.image_info import ImageInfoSet, ImageInfo

import warnings
from astropy.io.fits.verify import VerifyWarning
warnings.filterwarnings('ignore', category=VerifyWarning)

class ATLASParser():
    
    def __init__(self, working_path, image_type='science', crop_params=None):
        
        # paths
        self.working_path = working_path
        if not os.path.exists(self.working_path): os.makedirs(self.working_path)
        if not os.path.exists(os.path.join(self.working_path, 'ori/')): os.makedirs(os.path.join(self.working_path, 'ori/'))
        self.ori_image_path = os.path.join(self.working_path, 'ori/')
        if not os.path.exists(os.path.join(self.working_path, 'prep_images/')): os.makedirs(os.path.join(self.working_path, 'prep_images/'))
        self.kbmod_image_path = os.path.join(self.working_path, 'prep_images/')

        # config parameters
        self.rebinning_factor = 8
        self.image_type = image_type

        if crop_params is not None:
            self.crop = True
            self.crop_params = crop_params
        else:
            self.crop = False

    def oversample_array(self, array):
        rows = np.repeat(array, self.rebinning_factor, axis=0)
        return np.repeat(rows, self.rebinning_factor, axis=1)      
      
    def correct_header(self, header):
        header['CTYPE1'] = 'RA---TAN-SIP'
        header['CTYPE2'] = 'DEC--TAN-SIP'
        for key in list(header.keys()):
            if key.startswith('PV'): header.remove(key)
        header['RADESYSa'] = header['RADECSYS']
        del header['RADECSYS']
        header['OBS-LAT'] = header['SITELAT']
        header['OBS-LONG'] = header['SITELONG']
        header['OBS-ELEV'] = header['SITEELEV']
        if header['OBSID'][:2] == '01': header['OBSERVAT'] = 'T05'
        elif header['OBSID'][:2] == '02': header['OBSERVAT'] = 'T08'
        elif header['OBSID'][:2] == '03': header['OBSERVAT'] = 'W68'
        elif header['OBSID'][:2] == '04': header['OBSERVAT'] = 'M22'
        header['DATE-AVG'] = (datetime.strptime(header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S%z') +
                              timedelta(seconds=header['EXPTIME']/2)).strftime('%Y-%m-%dT%H:%M:%S')
        header['MJD-AVG'] = header['MJD-OBS'] + header['EXPTIME']/2/86400
        header["IDNUM"] = header['OBSID']
        header['PIXSCALE'] = header['RP_ASCL']
        return header
    
    def get_image(self, file, save=False):

        with fits.open(file) as hdul:
            img = hdul[1].data
            mask = img#hdul[2].data
            var = img#self.oversample_array(hdul[3].data).astype(np.uint16)
            header = hdul[1].header#self.correct_header(hdul[1].header)
        lab_img = os.path.basename(file).replace('fits.fz','sci.fits')

        if self.image_type == 'diff':
            sci_file = file.replace('fits.fz', 'diff.fz')
            img = fits.getdata(sci_file,1)
            lab_img = os.path.basename(file).replace('fits.fz','diff.fits')

        if self.crop:
            center_coords = self.crop_params['center_coords']
            crop_size = self.crop_params['crop_size']
            img, _ = self.crop_image(img, header, center_coords, crop_size)
            mask,_ = self.crop_image(mask, header, center_coords, crop_size)
            var, header = self.crop_image(var, header, center_coords, crop_size)
            header = self.correct_header(header)


        if save:
            image = [fits.PrimaryHDU(img, header), fits.PrimaryHDU(mask, header), fits.PrimaryHDU(var, header)]
            lab = [lab_img,
                   os.path.basename(file).replace('fits.fz','mask.fits'),
                   os.path.basename(file).replace('fits.fz','var.fits')]
            for h in range(len(image)):
                image[h].writeto(os.path.join(self.kbmod_image_path, lab[h]), overwrite=True, output_verify='silentfix')
        
        return img, header, mask, var
    
    def load_images(self, im_list):

        print("---------------------------------------")
        print("Loading Images")
        print("---------------------------------------")

        img_info = ImageInfoSet()
        images = []
        visit_times = []

        for im in im_list:

            sci, header, mask, var = self.get_image(im)

            header_info = ImageInfo()
            header_info.populate_from_header(header)

            time_stamp = header['MJD-AVG']
            psf = kb.psf(header['FWHM'])

            img = kb.layered_image(kb.raw_image(sci), 
                                    kb.raw_image(var),
                                    kb.raw_image(mask),
                                    time_stamp, psf)
            images.append(img)
            visit_times.append(time_stamp)
            img_info.append(header_info)
            

        stack = kb.image_stack(images)

        # Create a list of visit times and visit times shifted to 0.0.
        img_info.set_times_mjd(np.array(visit_times))
        times = img_info.get_zero_shifted_times()
        stack.set_times(times)
        print("Times set", flush=True)

        return stack, img_info
    '''

    def load_images_std(self, im_list):

        print("---------------------------------------")
        print("Loading Images")
        print("---------------------------------------")

        img_info = ImageInfoSet()
        images = []
        visit_times = []
        
        nnew = 0
        stdnew = 0
        avgnew = 0

        SCI = []; VAR = []; MASK = []
        TS = []; PSF = []

        for i,im in enumerate(im_list):

            sci, header, mask, var = self.get_image(im)

            header_info = ImageInfo()
            header_info.populate_from_header(header)

            time_stamp = header['MJD-AVG']
            psf = kb.psf(header['FWHM'])

            SCI.append(sci)
            VAR.append(var)
            MASK.append(mask)
            TS.append(time_stamp)
            PSF.append(psf)

            # img = kb.layered_image(kb.raw_image(sci), 
            #                         kb.raw_image(var),
            #                         kb.raw_image(mask),
            #                         time_stamp, psf)
            # images.append(img)
            visit_times.append(time_stamp)
            img_info.append(header_info)

            nnew, avgnew, stdnew = self.dyn_avgstd(sci, nnew, avgnew, stdnew)

        stdnew[stdnew==0] = 1e5

        for i in range(len(SCI)):
            img = kb.layered_image(kb.raw_image(SCI[i]/stdnew),
                                    kb.raw_image(VAR[i]/stdnew),
                                    kb.raw_image(MASK[i]),
                                    TS[i], PSF[i])
            images.append(img)


        stack = kb.image_stack(images)

        # Create a list of visit times and visit times shifted to 0.0.
        img_info.set_times_mjd(np.array(visit_times))
        times = img_info.get_zero_shifted_times()
        stack.set_times(times)
        print("Times set", flush=True)

        return stack, img_info
    '''
    
    @staticmethod
    def crop_image(im, header, coords, size=1000):
        wcs = WCS(header)
        center_coord = SkyCoord(coords[0], coords[1], unit=(u.deg, u.deg), frame='icrs')
        center_px = wcs.all_world2pix(center_coord.ra, center_coord.dec, 0)
        cutout = Cutout2D(im, center_px, size*u.pixel, wcs=wcs)
        cropped_data = cutout.data
        cropped_wcs = cutout.wcs
        header.update(cropped_wcs.to_header())
        return cropped_data, header
    

    @staticmethod
    def dyn_avgstd(valuenew, nold, avgold, stdold):
        nnew = nold + (valuenew != 0).astype(np.int)
        if np.sum(nold) == 0:
            avgnew = np.asarray(valuenew, dtype=np.double)
            stdnew = np.zeros_like(valuenew, dtype=np.double)
        else:
            avgnew = avgold + (valuenew - avgold) / nnew
            stdnew = np.sqrt(
                nold/nnew * stdold**2 +
                (valuenew - avgnew) * (valuenew - avgold) / nnew
                )
        return nnew, avgnew, stdnew