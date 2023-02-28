import os
import math
import time
import json
import numpy as np
import pandas as pd
import multiprocess as mp
# mp.set_start_method("spawn")

from .chromsizes import chromsizes


# convert a hilbert position from one order to another
def hilbert_pos_to_order(pos, fro, to):
    return math.floor(pos / pow(4, fro - to))

class HilbertGenome:
    def __init__(self, 
                 assembly='hg38',
                 output_dir='./data',
                 data=None,
                 name="order",
                 aggregator=lambda x: x,    #aggregator should always return a scalar or a numpy array
                 accessor=lambda x: x,
                 dtype="int8",
                 signal_resolution=1,
                 missing_value=-1, 
                ):
        self.assembly = assembly
        self.genome = chromsizes[assembly]
        self.max_order = 16 # this fits 3 billion bases in 4 billion points

        # all of the data we write out will 
        self.output_dir = output_dir
        if not os.path.exists(output_dir): 
            os.makedirs(output_dir)
        # this is a somewhat temporary directory 
        # where the files for each chromosome will be generated
        self.bytes_dir = os.path.join(output_dir, 'chromosome_bytes')
        if not os.path.exists(self.bytes_dir): 
            os.makedirs(self.bytes_dir)

        # we expect the data to be a dictionary of chromosomes
        # where each chromosome entry is an iterable of BED compatible entries
        self.data = data
        # actually we expect the data to be a pandas dataframe


        self.accessor = accessor
        self.aggregator = aggregator
        # TODO: lets get the width of the aggregated data with a subset here
        self.width = self.aggregator(data).size
        self.dtype = dtype
        self.name = name

        self.signal_resolution = signal_resolution
        self.missing_value = missing_value
        # TODO: print out a bunch more helpful upfront info
        print("data width", self.width)
        for order in range(6, self.max_order):
            individual = self.get_hilbert_genome_resolution(order) <= self.signal_resolution
            print("order", order, "individual", individual)

    def generate_order(self, order):
        """ Generate all the files for an order in parallel then combine them"""    
        start_time = time.time()
        files = self.files_for_order(order)
        print("generating", len(files), "files")
        for file in files:
            self.aggregate_parallel(file)
        print("order", order, "done in: %s" % (time.time() - start_time))
        print("combining files")
        self.combine_bytes_files(files, order)
        self.combine_meta_files(files, order)
        print("done with", order)


    def aggregate_parallel(self, params):
        start_time = time.time()
        print(params['chromosome'], params['order'], params['hstart'], "starting")
        v = self.aggregate_range_region(**params)
        self.write_aggregate_range_region(v, **params)
        print(params['chromosome'], params['order'], params['hstart'], "done in: %s" % (time.time() - start_time))

    def aggregate_chromosome(self, chromosome, order, chunk_size=4**11):
        genome = self.genome
        gstart = genome.sizes_acc[chromosome]
        gstop = genome.sizes_acc[chromosome] + genome.sizes[chromosome]
        hbstart = self.get_hilbert_from_genome(gstart)
        hbstop = self.get_hilbert_from_genome(gstop)
        
        hstart = hilbert_pos_to_order(hbstart, 16, order)
        hstop = hilbert_pos_to_order(hbstop, 16, order)
        hlen = hstop - hstart

        params = []
        if hlen > chunk_size:
            # break up our hlen into chunks
            nchunks = math.ceil(hlen/chunk_size)
            # print("chunks", nchunks)
            for i in range(0, nchunks):
                cstart = hstart + chunk_size * i
                cstop = cstart + chunk_size
                if cstop > hstop:
                    cstop = hstop
                params.append({
                    'chromosome': chromosome,
                    'order': order,
                    'hstart': cstart,
                    'hstop': cstop,
                    'gstart': gstart, # this stays the same
                })
                
        else:
            params.append({
                'chromosome': chromosome,
                'order': order,
                'hstart': hstart,
                'hstop': hstop,
                'gstart': gstart,
            })
        return params

    # for a given range of hilbert indices, aggregate the data
    def aggregate_range_region(self, chromosome, order, hstart, hstop, gstart):
        dtype = self.dtype
        # lets trim down the data to just the samples we are interested in
        df = self.data
        start = self.get_genome_from_hilbert_order(hstart, order) - gstart
        stop = self.get_genome_from_hilbert_order_reverse(hstop, order) - gstart
        # we will filter down to all the samples who start within the range of our hilbert range
        samples = df[(df['chromosome'] == chromosome) & (df['start'] >= start) & (df['start'] <= stop)]

        # if there are more hilbert indices than samples we dont need to aggregate
        individual = self.get_hilbert_genome_resolution(order) <= self.signal_resolution

        # the dimensions of our resulting array are determined by
        # the stretch of the hilbert indices we are working with
        # and the width determined by the aggregator
        rows = hstop - hstart
        width = self.width
        v = np.full((rows, width), self.missing_value, dtype=np.dtype(dtype))
        if(len(samples) == 0):
            return v

        if individual:
            # lets loop thru the samples and assign them to the hilbert indices
            for idx,sample in samples.iterrows():
                # get the hilbert index for the start of this sample
                hbstart = self.get_hilbert_from_genome(sample['start'] + gstart)
                # convert to the order we are working with
                shstart = hilbert_pos_to_order(hbstart, 16, order)
                # get the hilbert index for the end of this sample
                hbstop = self.get_hilbert_from_genome(sample['end'] + gstart)
                # convert to the order we are working with
                shstop = hilbert_pos_to_order(hbstop, 16, order)

                istart = shstart - hstart
                istop = shstop - hstop
                if(istart < 0):
                    istart = 0
                if(istop > rows):
                    istop = rows
                v[istart:istop] = self.accessor(sample)
        else:
            # we copy the dataframe and add a column for the hilbert index
            samples_copy = samples.copy()
            samples_copy["hpos"] = samples_copy.apply(lambda x: 
                self.get_hilbert_from_genome_order(x['start'] + gstart, order)
            , axis=1)

            # then we group by the hilbert index and aggregate
            result = samples_copy.groupby("hpos").apply(self.aggregator)
            # if width is 1 there will always only be a single value
            for pos in result.index:
                i = pos - hstart
                if(i < 0 or i >= rows):
                    continue
                v[i] = result.loc[pos]

        return v

    def write_aggregate_range_region(self, v, chromosome, order, hstart, hstop, gstart):
        dtype = self.dtype
        fn = self.get_bytes_filename(chromosome, order, hstart, dtype)
        with open(fn, 'wb') as f:
            f.write(v.tobytes()) 
        
        # lets calculate the minimum and maximum values
        # and save them in a metadata json file
        min = np.amin(v, axis=0).tolist()
        max = np.amax(v, axis=0).tolist()

        meta = {'min': min, 'max': max, 'dtype': dtype, 'width': self.width }
        mfn = self.get_meta_filename(chromosome, order, hstart)
        with open(mfn, 'w') as f:
            json.dump(meta, f)


    def files_for_order(self, order):
        files = []
        for c in self.genome.chromosomes:
            files += self.aggregate_chromosome(c, order)
        return files

    def get_bytes_filename(self, chromosome, order, hstart, dtype):
        return os.path.join(self.bytes_dir, self.name + "_" + chromosome + "_" + str(order) + "_" + str(hstart) + "." + dtype)

    def get_meta_filename(self, chromosome, order, hstart):
        return os.path.join(self.bytes_dir, self.name + "_" + chromosome + "_" + str(order) + "_" + str(hstart) + ".json")


    def combine_bytes_files(self, files, order):
        paths = [self.get_bytes_filename(f['chromosome'], f['order'], f['hstart'], self.dtype) for f in files]
        print("Combining files", paths)
        out = os.path.join(self.output_dir, self.name + "_" + str(order) + "." + self.dtype)
        with open(out, "wb") as of:
            for p in paths:
                with open(p, "rb") as rf:
                    of.write(rf.read())
        print("all done writing", out)

    def combine_meta_files(self, files, order):
        paths = [self.get_meta_filename(f['chromosome'], f['order'], f['hstart']) for f in files]
        print("Combining files", paths)
        out = os.path.join(self.output_dir, self.name + "_" + str(order) + ".json")
        with open(out, "w") as of:
            mins = []
            maxs = []
            dtype = None
            for p in paths:
                with open(p, "r") as rf:
                    m = json.loads(rf.read())
                    dtype = m['dtype'] # this will always be the same
                    width = m['width'] # this will always be the same
                    mins.append(np.array(m['min']))
                    maxs.append(np.array(m['max']))
            
            meta = {
                'min': np.array(mins).min(axis=0).tolist(),
                'max': np.array(maxs).max(axis=0).tolist(),
                'dtype': dtype,
                'width': width
            }
            of.write(json.dumps(meta))
        print("all done!")



    # spread our genome positions evenly across the 16th order hilbert curve
    def spread_scale(self, gpos):
        i_min = 0
        i_max = self.genome.total_size
        o_min = 0
        o_max = pow(4, self.max_order)
        return round(o_min + (gpos - i_min) * ((o_max - o_min) / (i_max - i_min)))

    # get the genome position for a given hilbert position (may return duplicates)
    def spread_scale_invert(self, pos):
        i_min = 0
        i_max = self.genome.total_size
        o_min = 0
        o_max = pow(4, self.max_order)
        return round((pos - o_min) / ((o_max - o_min) / (i_max - i_min)) + i_min)

    # get the global genome position from something like "chr2:12345"
    def get_genome_from_coord(self, start):
        chr = start.split(":")[0]
        offset = int(start.split(":")[1])
        return self.genome.sizes_acc[chr] + offset

    def get_hilbert_from_genome(self, gpos):
        return self.spread_scale(gpos)

    def get_hilbert_from_genome_order(self, gpos, order):
        hbp = self.get_hilbert_from_genome(gpos)
        return hilbert_pos_to_order(hbp, self.max_order, order)


    def get_genome_from_hilbert_base(self, pos):
        gpos = self.spread_scale_invert(pos)
        if pos != self.spread_scale(gpos): 
            return None
        return gpos

    # get the first base at the maximum order resolution 
    # that comes after the starting pos and before the next pos at the current order
    def get_genome_from_hilbert_order(self, pos, order):
        start = hilbert_pos_to_order(pos, order, self.max_order)
        found = self.get_genome_from_hilbert_base(start)
        # need to scan everything at the highest order
        # but only until we find something which should usually happen in a few loops
        if order == 16:
            # if we are at the max order, we actually want to search a little further
            order = 15
        for l in range(pow(4, self.max_order - order)):
            if found is not None:
                return found
            found = self.get_genome_from_hilbert_base(start+l)

    # get the last base at the maximum order resolution 
    # that comes after the starting pos and before the next pos at the current order
    def get_genome_from_hilbert_order_reverse(self, pos, order):
        start = hilbert_pos_to_order(pos+1, order, self.max_order) - 1
        found = self.get_genome_from_hilbert_base(start)
        # need to scan everything at the highest order
        # but only until we find something which should usually happen in a few loops
        if order == 16:
            # if we are at the max order, we actually want to search a little further
            order = 15
        for l in range(pow(4, self.max_order - order)):
            if found is not None:
                return found
            found = self.get_genome_from_hilbert_base(start-l)

    # get the approximate number of base pairs per hilbert celll at a given order
    def get_hilbert_genome_resolution(self, order):
        return self.get_genome_from_hilbert_order_reverse(0, order) + 1