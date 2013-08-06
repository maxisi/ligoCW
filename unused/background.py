import tables as tb
import numpy as np
from analysis import heterodyne

def rhsave(data, freq, det):
    
    rhB_dir = 'files/background/' + det + '.h5'
    
    with tb.open_file(rhB_dir, mode = 'w', title = det + 'background') as h5file:

        class Data(tb.IsDescription):
            rh_t = tb.Int32Col()        # time
            rh_r = tb.Float64Col()      # real part
            rh_i = tb.Float64Col()      # imag part

    
#         print 'created file'
    
        for psr in data.columns:
    
            psrname = psr.replace('+','p').replace('-','m')
        
            group = h5file.create_group('/', psrname, 'PSR' + psrname)
        
#             print 'created group %s' % psrname
        
            rh = heterodyne.het(data[psr], freq)
            
#             print 'heterodyned'
            
            for f_i in range(0, len(freq)):
#                 print 'creating table...',
                
                table_name = 'f' + str(f_i)
                
                table = h5file.create_table(group, table_name, Data, "Inst" + str(f_i) )
                
                rhf = rh[table_name]
                
#                 print table_name
                
                entry = table.row
                
                for t in rhf.index:
            
                    entry['rh_t'] = t
                    entry['rh_r'] = rhf[t].real
                    entry['rh_i'] = rhf[t].imag
                
                    entry.append()
    
                table.flush()