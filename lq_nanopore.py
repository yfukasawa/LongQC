import os, sys, time, h5py, json
import shutil, random, tarfile, logging
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy             as np

from multiprocessing import Pool, Manager, Lock
from operator        import itemgetter

def get_flowcell_coord():
    # for the case we have multiple layout
    return _c2cor_r94_r95()

def _cor2c_r94_r95():
    layout = np.zeros((32,16), dtype=int)

    asc  = [33, 481, 417, 353, 289, 225, 161, 97]
    desc = [1, 449, 385, 321, 257, 193, 129, 65]
    for i,num in enumerate(asc):
        for j in np.arange(4):
            for z, c in enumerate(np.arange(num+8*j, num+8*j+8)):
                layout[i*4+j][z] = c
    for i,num in enumerate(desc):
        for j in np.arange(4):
            for z, c in enumerate(np.arange(num+8*j, num+8*j+8)):
                layout[i*4+j][15-z] = c

    return layout

def _c2cor_r94_r95():
    layout = [0] * 513

    asc  = [33, 481, 417, 353, 289, 225, 161, 97]
    desc = [1, 449, 385, 321, 257, 193, 129, 65]

    for i,num in enumerate(asc):
        for j in np.arange(4):
            for z, c in enumerate(np.arange(num+8*j, num+8*j+8)):
                layout[c] = (i*4+j,z)
    for i,num in enumerate(desc):
        for j in np.arange(4):
            for z, c in enumerate(np.arange(num+8*j, num+8*j+8)):
                layout[c] = (i*4+j,15-z)
    layout[0] = None

    return layout

def list_fast5_targz(d):
    if not os.path.isdir(d):
        return None

    list_ftg = []
    list = os.listdir(d)

    for i in list:
        p = os.path.join(d, i)
        if os.path.isdir(p):
            continue

        if p.endswith("tar.gz"):
            list_ftg.append(p)

    return list_ftg

def get_members_from_tar(fname):
    tar = tarfile.open(fname)
    return tar.getmembers()

#h5py cannot accept file handle, so we have to get actual files first.
def extract_tar(fname, base_path=''):
    #d   = os.path.dirname(os.path.abspath(fname))
    tar = tarfile.open(fname)
    tar.extractall(base_path)
    tar.close()

#https://community.nanoporetech.com/posts/pulling-time-to-each-read?search_term=start_time
def list_fast5_files(d, logger):
    if not os.path.isdir(d):
        logger.info("%s is not a dir" % d)
        return None

    list_fast5 = []
    list = os.listdir(d)

    for i in list:
        p = os.path.join(d, i)
        if os.path.isdir(p):
            logger.info("Go into subdir for listing fast5: %s" % p)
            _list = os.listdir(p)
            for j in _list:
                _p = os.path.join(p, j)
                
                if _p.endswith("fast5"):
                    list_fast5.append(_p)

        if p.endswith("fast5"):
            list_fast5.append(p)

    return list_fast5

def open_fast5(path):
    from lq_utils import eprint
    try:
        f = h5py.File(path, 'r')
    except OSError as e:
        eprint(e)
        return None
    return f

def list_toplevel(f):
    return list(f.keys())

def get_fastq_from_multi_fast5(f, rn):
    # fastq binary string is stored under /Analyses/Basecall_1D_000/BaseCalled_template/Fastq
    # empty tuple index means 'scalar' access in a 'dataset' for h5py
    return f[rn]['Analyses']['Basecall_1D_000']['BaseCalled_template']['Fastq'][()].decode('ascii')

def get_channel_id(f):
    return int(f['/UniqueGlobalKey']['channel_id'].attrs['channel_number'])

def get_sampling_rate(f):
    return int(f['/UniqueGlobalKey']['channel_id'].attrs['sampling_rate'])

def get_read_nodename(f):
    keys = f['Raw/Reads'].keys()

    if len(keys) == 0:
        # error
        pass

    return list(keys)[0]

def get_flowcell(f):
    return f['/UniqueGlobalKey']['context_tags'].attrs['flowcell_type']

def get_kit(f):
    # sequencing_kit
    return f['/UniqueGlobalKey']['context_tags'].attrs['sequencing_kit']

def get_start_time(f):
    n   = get_read_nodename(f)
    s_t = f['Raw/Reads'][n].attrs['start_time']
    return int(s_t/get_sampling_rate(f))

def get_duration(f):
    n = get_read_nodename(f)
    d = f['Raw/Reads'][n].attrs['duration']
    return int(d/get_sampling_rate(f))

def wrapper(f5):
    fast5 = open_fast5(f5)
    if fast5 == None:
        return None
    c_id  = get_channel_id(fast5) -1 
    s_t   = get_start_time(fast5)
    durat = get_duration(fast5)
    fc  = get_flowcell(fast5)
    kit = get_kit(fast5)
    fast5.close()

    return ( c_id, (s_t, s_t+durat), fc, kit )

def process_uncompressed_single(l, logger, n_channel=512):
    logger.info("File list loaded. Size = %d" % len(l))

    # single core (for test)
    start_time             = time.time()
    channel_wise_occupancy = [None] * n_channel
    occ                    = []

    for f5 in l:
        fast5 = open_fast5(f5)
        if fast5 == None:
            continue
        c_id = get_channel_id(fast5) -1 
        if channel_wise_occupancy[c_id] == None:
            channel_wise_occupancy[c_id] = []

        s_t = get_start_time(fast5)
        channel_wise_occupancy[c_id].append((s_t, s_t+get_duration(fast5)))

        fast5.close()

def process_uncompressed_multi(p, l, bag, fcs, kits, logger):
    logger.info("File list loaded. Size = %d" % len(l))

    #bag = [set() for _ in range(512)]

    rtn_state = p.map(wrapper, l)

    for t in rtn_state:
        if t == None:
            continue
        bag[t[0]].add(t[1])
        fcs.add(t[2])
        kits.add(t[3])

        """
        # assertion
        if len(bag) != len(channel_wise_occupancy):
            print("Error")
        
        cnt = 0
        for i in range(len(bag)):
            if len(bag[i]) == 0 and channel_wise_occupancy[i] == None:
                continue

            if bag[i] != set(channel_wise_occupancy[i]):
                print(i)
                print(bag[i], " ", channel_wise_occupancy[i])
                cnt +=1

        print(cnt)
        """

def run_platformqc(platform, data_path, output_path, *, suffix=None, n_channel = 512, n_process=15):
    THRESHOLD_INACTIVE = 0.0025

    ld = os.path.join(output_path, "log")
    pd = os.path.join(output_path, "fig")
    if not os.path.isdir(ld):
        os.makedirs(ld, exist_ok=True)

    if not os.path.isdir(pd):
        os.makedirs(pd, exist_ok=True)
    if not suffix:
        suffix = ""
    else:
        suffix = "_" + suffix
    log_path  = os.path.join(ld, "log_ont_platform" +  suffix + ".txt")
    plot_path = os.path.join(pd, "fig_ont_platform" +  suffix + ".png")
    json_path = os.path.join(output_path, "QC_vals_" + platform + suffix + ".json")
    # json
    tobe_json = {}

    ### logging conf ###
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    fh = logging.FileHandler(log_path, 'w')
    sh = logging.StreamHandler()

    formatter = logging.Formatter('%(module)s:%(asctime)s:%(lineno)d:%(levelname)s:%(message)s')
    fh.setFormatter(formatter)
    sh.setFormatter(formatter)

    logger.addHandler(sh)
    logger.addHandler(fh)
    #####################

    logger.info("Started %s platform QC for %s" % (platform, data_path))
    logger.info("The num of channel:%d, the num of process:%d" % (n_channel, n_process))

    l    = list_fast5_files(data_path, logger)
    ltgz = list_fast5_targz(data_path)

    if len(l) == 0 and len(ltgz) == 0:
        logger.warning("No fast5 or compressed file in the given path: %s" % data_path)
        return 1

    if len(l) > 0 and len(ltgz) > 0:
        logger.warning("The given path is a mixture of compressed and uncompressed files. Check: %s" % data_path)
        return 1

    bag  = [set() for _ in range(n_channel)]
    fcs  = set()
    kits = set()
    occ  = []

    p = Pool(processes=n_process)

    if len(l) == 0:
        logger.info("There is no fast5 file in the given dir: %s" % data_path)
        #print(l)
        for f in ltgz:
            logger.info("Process a compressed file: %s" % f)
            base_dir = os.path.dirname( os.path.abspath(f) )
            sub_dir  = os.path.basename(f).replace(".tar.gz", '')
            extract_tar(f, base_path=base_dir) # this method extract all members at the same dir of f
            _l = list_fast5_files( os.path.join(base_dir, sub_dir), logger )
            process_uncompressed_multi(p, _l, bag, fcs, kits, logger)
            shutil.rmtree( os.path.join(base_dir, sub_dir) )
    else:
        process_uncompressed_multi(p, l, bag, fcs, kits, logger)

    p.close()
    logger.info("All of fast5 files were loaded.")

    tobe_json['Sequencing kit'] = str(", ".join(str(s.decode("utf-8")) for s in kits))
    tobe_json['Flowcell'] = str(", ".join(str(s.decode("utf-8")) for s in fcs))

    logger.info("Aggregating pore running time.")
    max = -1
    channel_wise_cnt = [0] * n_channel
    for i in range(len(bag)):
        s = sorted(bag[i], key=itemgetter(0,1))
        bag[i] = s
        if len(s) > 0 and s[-1][1] > max:
            max = s[-1][1]

    for i in range(1, max+1):
        cnt = 0
        for j in range(len(bag)):
            if len(bag[j]) < 1:
                continue

            if bag[j][0][0] <= i and i <= bag[j][0][1]:
                cnt += 1
                channel_wise_cnt[j] += 1

            elif bag[j][0][1] < i:
                bag[j].pop(0)

        occ.append(cnt/len(bag))

    logger.info("Aggregation finished.")

    tobe_json['Sequencing time in seconds'] = int(max)
    tobe_json['The time reached maximum active pore rate'] = int(np.argmax(occ))
    tobe_json['The maximum active pore rate'] = float(np.max(occ))

    for i in range(n_channel):
        channel_wise_cnt[i] /= max

    tobe_json['The fraction of inactive pores'] = float(np.where(np.array(channel_wise_cnt) < THRESHOLD_INACTIVE)[0].shape[0]/n_channel)

    y = np.arange(0, 33)
    x = np.arange(0, 17)
    X,Y = np.meshgrid(x,y)
    Z = np.zeros((33,17), dtype=float)
    c2cor = get_flowcell_coord()

    logger.info("Generating plot 1.")

    for (c, cor)  in enumerate(c2cor):
        if cor is None:
            continue
        Z[cor[0]][cor[1]] = channel_wise_cnt[c-1]

    #plt.figure(figsize=(5,4))
    plt.subplot(3,1,1)
    plt.plot(occ)
    plt.grid(True)

    plt.xlabel('Elapsed time in seconds')
    plt.ylabel('Active channel rate')

    for i in np.arange(1, max+1, 28800): # 8 hours
        if i == 1:
            continue
        plt.axvline(x=i, linestyle='dashed', linewidth=1, color='blue', alpha=0.8)

    logger.info("Generating plot 2.")
    plt.subplot(3,1,2)
    plt.pcolor(X, Y, Z, cmap='RdBu')
    plt.colorbar()
    plt.tight_layout()
    plt.title("Pore activity mapped on the actual layout")
    plt.contour(X, Y, Z, levels=[THRESHOLD_INACTIVE], linewidths=2, linestyles='dashed')
    plt.pink()

    logger.info("Generating plot 3.")
    plt.subplot(3,1,3)
    plt.hist(channel_wise_cnt, color='blue', bins=100)
    plt.xlabel('Channel wise activity rate')
    plt.ylabel('Frequency')

    plt.subplots_adjust(hspace=1.0)

    plt.savefig(plot_path, bbox_inches="tight")
    #plt.savefig(plot_path)
    #plt.show()
    logger.info("Plots were saved to %s." % plot_path)

    with open(json_path, "w") as f:
        logger.info("Quality measurements were written into a JSON file: %s" % json_path)
        json.dump(tobe_json, f, indent=4)

# test
if __name__ == "__main__":
    d = "/home/fukasay/rawdata/ont/20171029_1450_20171029_ecoli_1D_square_test/fast5/"
    run_platformqc(d, "/home/fukasay/analyses/longQC/ont_platform_test_10291d/", n_process=5)
