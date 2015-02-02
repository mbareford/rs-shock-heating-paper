import array
import struct
import os


def get_index(c_val, c_min, c_max, cn):
  
    c_i = 0
    
    if (c_val > c_max):
        c_i = cn - 1
    elif (c_val >= c_min):
        dc = float(abs(c_max-c_min))/float(cn-1)
        c = c_min
                        
        while (c_i < cn):
        
            if (c > c_val):
                if (abs(c_val-c) > abs(c_val-(c-dc))):    
                    c_i = c_i - 1
                break
                    
            c_i = c_i + 1
            c = c_min + (dc*c_i)
                  
    return c_i  
  
    
def get_field_buffer(val, val_format, val_size):    
        
    field = array.array('b',[0])*val_size
    struct.pack_into(val_format, field, 0, val)
    
    return field


def get_string_field_buffer(val, val_len):    
        
    field = array.array('b',[0])*val_len
    val_format = str(eval('val_len')) + 's'
    struct.pack_into(val_format, field, 0, val)
    
    return field


def get_array_field_buffer(val, val_format):    
    
    return struct.pack(val_format % len(val), *val)
                

# python -c 'import lare3d_cfd as cfd; cfd.slice("../../data/256x256x512/0095.cfd.org", "../../data/256x256x512/0095.cfd", -1.0, 1.0)'  
def slice(cfd_path,cfd_seg_path,z_min,z_max):
                           
    char_size = 1
    int_size = 4
    long_size = 8
    dbl_size = 8
    hdr_tag_len = 3    
    

    cfd = open(cfd_path, 'rb')
    cfd_seg = open(cfd_seg_path, 'wb')
    
    
    # process file header
    #################################################################
    hdr_tag = cfd.read(hdr_tag_len)        
    hdr_size = struct.unpack_from('i', cfd.read(int_size))[0]                
    blk_hdr_size = struct.unpack_from('i', cfd.read(int_size))[0]
    cfd_ver = struct.unpack_from('i', cfd.read(int_size))[0]
    cfd_rev = struct.unpack_from('i', cfd.read(int_size))[0]
    max_str_len = struct.unpack_from('i', cfd.read(int_size))[0]
    nblocks = struct.unpack_from('i', cfd.read(int_size))[0]
    
    
    print 'Header tag: ', hdr_tag
    print 'Header size: ', hdr_size
    print 'Block header size: ', blk_hdr_size
    print 'CFD version: ', cfd_ver
    print 'CFD revision: ', cfd_rev
    print 'Max str len: ', max_str_len
    print 'nblocks: ', nblocks    
    
    
    cfd_seg.write(get_string_field_buffer(hdr_tag,hdr_tag_len))
    cfd_seg.write(get_field_buffer(hdr_size,'i',int_size))
    cfd_seg.write(get_field_buffer(blk_hdr_size,'i',int_size))
    cfd_seg.write(get_field_buffer(cfd_ver,'i',int_size))
    cfd_seg.write(get_field_buffer(cfd_rev,'i',int_size))
    cfd_seg.write(get_field_buffer(max_str_len,'i',int_size))
    cfd_seg.write(get_field_buffer(nblocks,'i',int_size))        
    #################################################################

    
    print ' '
    
    
    # process first block (time)
    ######################################################################
    tm_blk_name = cfd.read(max_str_len)
    tm_blk_class = cfd.read(max_str_len)
    tm_blk_type = struct.unpack_from('i', cfd.read(int_size))[0]
    tm_blk_len = struct.unpack_from('l', cfd.read(long_size))[0]
    tm_blk_meta_len = struct.unpack_from('l', cfd.read(long_size))[0]
       
    cycle = struct.unpack_from('i', cfd.read(int_size))[0]
    time = struct.unpack_from('d', cfd.read(dbl_size))[0]
    
    
    print 'Time block name: ', tm_blk_name
    print 'Time block class: ', tm_blk_class
    print 'Time block type: ', tm_blk_type
    print 'Time block len: ', tm_blk_len
    print 'Time block meta len: ', tm_blk_meta_len
    cfd_seg_path
    print 'Cycle: ', cycle
    print 'Time: ', time
    
    
    cfd_seg.write(get_string_field_buffer(tm_blk_name,max_str_len))
    cfd_seg.write(get_string_field_buffer(tm_blk_class,max_str_len))    
    cfd_seg.write(get_field_buffer(tm_blk_type,'i',int_size))
    cfd_seg.write(get_field_buffer(tm_blk_len,'l',long_size))
    cfd_seg.write(get_field_buffer(tm_blk_meta_len,'l',long_size))
    
    cfd_seg.write(get_field_buffer(cycle,'i',int_size))
    cfd_seg.write(get_field_buffer(time,'d',dbl_size))
    ######################################################################
    
    
    print ' '
    

    # process second block (grid)
    ##############################################################################################################
    grd_blk_name = cfd.read(max_str_len)
    grd_blk_class = cfd.read(max_str_len)
    grd_blk_type = struct.unpack_from('i', cfd.read(int_size))[0]
    grd_blk_len = struct.unpack_from('l', cfd.read(long_size))[0]
    grd_blk_meta_len = struct.unpack_from('l', cfd.read(long_size))[0]
    
    grd_type = struct.unpack_from('i', cfd.read(int_size))[0]
    grd_dim = struct.unpack_from('i', cfd.read(int_size))[0]
    grd_field_size = struct.unpack_from('i', cfd.read(int_size))[0]    
    grd_nx = struct.unpack_from('i', cfd.read(int_size))[0]
    grd_ny = struct.unpack_from('i', cfd.read(int_size))[0]
    grd_nz = struct.unpack_from('i', cfd.read(int_size))[0]    
    grd_x_min = struct.unpack_from('d', cfd.read(dbl_size))[0]
    grd_x_max = struct.unpack_from('d', cfd.read(dbl_size))[0]
    grd_y_min = struct.unpack_from('d', cfd.read(dbl_size))[0]
    grd_y_max = struct.unpack_from('d', cfd.read(dbl_size))[0]
    grd_z_min = struct.unpack_from('d', cfd.read(dbl_size))[0]
    grd_z_max = struct.unpack_from('d', cfd.read(dbl_size))[0]
    
    grd_x = struct.unpack_from(str(eval('grd_nx')) + 'd', cfd.read(dbl_size*grd_nx))
    grd_y = struct.unpack_from(str(eval('grd_ny')) + 'd', cfd.read(dbl_size*grd_ny))
    
    
    if (z_min == grd_z_min and z_max == grd_z_max):
        grd_z = struct.unpack_from(str(eval('grd_nz')) + 'd', cfd.read(dbl_size*grd_nz))
    else:
        z_min_i = get_index(z_min, grd_z_min, grd_z_max, grd_nz)
        z_max_i = get_index(z_max, grd_z_min, grd_z_max, grd_nz)
        cfd.seek(z_min_i*grd_field_size, 1)
        grd_z = struct.unpack_from(str(eval('z_max_i-z_min_i+1')) + 'd', cfd.read(grd_field_size*(z_max_i-z_min_i+1)))
        cfd.seek((grd_nz-z_max_i-1)*grd_field_size, 1)
        grd_nz = z_max_i-z_min_i+1
        grd_z_min = min(grd_z)
        grd_z_max = max(grd_z)
        grd_blk_meta_len = grd_blk_len + (grd_nx + grd_ny + grd_nz)*grd_field_size
    
    
    print 'Grid block name: ', grd_blk_name
    print 'Grid block class: ', grd_blk_class
    print 'Grid block type: ', grd_blk_type
    print 'Grid block len: ', grd_blk_len
    print 'Grid block meta len: ', grd_blk_meta_len
    
    print 'Grid type: ', grd_type
    print 'Grid dim: ', grd_dim
    print 'Grid field size: ', grd_field_size    
    print 'Grid nx: ', grd_nx
    print 'Grid ny: ', grd_ny
    print 'Grid nz: ', grd_nz    
    print 'Grid x min: ', grd_x_min
    print 'Grid x max: ', grd_x_max
    print 'Grid y min: ', grd_y_min
    print 'Grid y max: ', grd_y_max
    print 'Grid z min: ', grd_z_min
    print 'Grid z max: ', grd_z_max    
    
    #print 'Grid x: ', grd_x
    #print 'Grid y: ', grd_y
    #print 'Grid z: ', grd_z
    
    
    cfd_seg.write(get_string_field_buffer(grd_blk_name,max_str_len))
    cfd_seg.write(get_string_field_buffer(grd_blk_class,max_str_len))    
    cfd_seg.write(get_field_buffer(grd_blk_type,'i',int_size))
    cfd_seg.write(get_field_buffer(grd_blk_len,'l',long_size))
    cfd_seg.write(get_field_buffer(grd_blk_meta_len,'l',long_size))
    
    cfd_seg.write(get_field_buffer(grd_type,'i',int_size))
    cfd_seg.write(get_field_buffer(grd_dim,'i',int_size))
    cfd_seg.write(get_field_buffer(grd_field_size,'i',int_size))
    cfd_seg.write(get_field_buffer(grd_nx,'i',int_size))
    cfd_seg.write(get_field_buffer(grd_ny,'i',int_size))
    cfd_seg.write(get_field_buffer(grd_nz,'i',int_size))    
    cfd_seg.write(get_field_buffer(grd_x_min,'d',dbl_size))
    cfd_seg.write(get_field_buffer(grd_x_max,'d',dbl_size))
    cfd_seg.write(get_field_buffer(grd_y_min,'d',dbl_size))
    cfd_seg.write(get_field_buffer(grd_y_max,'d',dbl_size))
    cfd_seg.write(get_field_buffer(grd_z_min,'d',dbl_size))
    cfd_seg.write(get_field_buffer(grd_z_max,'d',dbl_size))
    
    cfd_seg.write(get_array_field_buffer(grd_x,'%sd'))
    cfd_seg.write(get_array_field_buffer(grd_y,'%sd'))
    cfd_seg.write(get_array_field_buffer(grd_z,'%sd'))
    ##############################################################################################################      
    
        
    for i in range(2, nblocks):
        print ' '
        
        # process block i
        #################################################################################################################################################
        var_blk_name = cfd.read(max_str_len)
        var_blk_class = cfd.read(max_str_len)
        var_blk_type = struct.unpack_from('i', cfd.read(int_size))[0]
        var_blk_len = struct.unpack_from('l', cfd.read(long_size))[0]
        var_blk_meta_len = struct.unpack_from('l', cfd.read(long_size))[0]
    
        var_type = struct.unpack_from('i', cfd.read(int_size))[0]
        var_dim = struct.unpack_from('i', cfd.read(int_size))[0]
        var_field_size = long(struct.unpack_from('i', cfd.read(int_size))[0])
    
        var_nx = long(struct.unpack_from('i', cfd.read(int_size))[0])
        var_ny = long(struct.unpack_from('i', cfd.read(int_size))[0])
        var_nz = long(struct.unpack_from('i', cfd.read(int_size))[0])        
        var_stag_x = struct.unpack_from('d', cfd.read(dbl_size))[0]
        var_stag_y = struct.unpack_from('d', cfd.read(dbl_size))[0]
        var_stag_z = struct.unpack_from('d', cfd.read(dbl_size))[0]    
        var_min = struct.unpack_from('d', cfd.read(dbl_size))[0]
        var_max = struct.unpack_from('d', cfd.read(dbl_size))[0]    
        var_mesh_name = cfd.read(max_str_len)
        var_mesh_class = cfd.read(max_str_len)
        
        
        if (z_min == grd_z_min and z_max == grd_z_max):
            var_data = struct.unpack_from(str(eval('var_nz*var_nx*var_ny')) + 'd', cfd.read(var_field_size*(var_nz*var_nx*var_ny)))
        else:                        
            cfd.seek(long(z_min_i)*var_nx*var_ny*var_field_size, 1)
            var_data = struct.unpack_from(str(eval('long(z_max_i-z_min_i+1)*var_nx*var_ny')) + 'd', cfd.read(var_field_size*(long(z_max_i-z_min_i+1)*var_nx*var_ny)))
            cfd.seek(long(var_nz-z_max_i-1)*var_nx*var_ny*var_field_size, 1)
            var_nz = z_max_i-z_min_i+1
            var_min = min(var_data)
            var_max = max(var_data)
            var_blk_meta_len = var_blk_len + (var_nx*var_ny*var_nz)*var_field_size
    
                    
        print 'Variable block name: ', var_blk_name
        print 'Variable block class: ', var_blk_class
        print 'Variable block type: ', var_blk_type
        print 'Variable block len: ', var_blk_len
        print 'Variable block meta len: ', var_blk_meta_len
        print 'Variable type: ', var_type
        print 'Variable dim: ', var_dim
        print 'Variable field size: ', var_field_size
        print 'Variable nx: ', var_nx
        print 'Variable ny: ', var_ny
        print 'Variable nz: ', var_nz    
        print 'Variable stagger x: ', var_stag_x
        print 'Variable stagger y: ', var_stag_y
        print 'Variable stagger z: ', var_stag_z
        print 'Variable min: ', var_min
        print 'Variable max: ', var_max
        print 'Variable mesh name: ', var_mesh_name
        print 'Variable mesh class: ', var_mesh_class
                
        
        cfd_seg.write(get_string_field_buffer(var_blk_name,max_str_len))
        cfd_seg.write(get_string_field_buffer(var_blk_class,max_str_len))    
        cfd_seg.write(get_field_buffer(var_blk_type,'i',int_size))
        cfd_seg.write(get_field_buffer(var_blk_len,'l',long_size))
        cfd_seg.write(get_field_buffer(var_blk_meta_len,'l',long_size))
                
        cfd_seg.write(get_field_buffer(var_type,'i',int_size))
        cfd_seg.write(get_field_buffer(var_dim,'i',int_size))
        cfd_seg.write(get_field_buffer(var_field_size,'i',int_size))
             
        cfd_seg.write(get_field_buffer(var_nx,'i',int_size))
        cfd_seg.write(get_field_buffer(var_ny,'i',int_size))
        cfd_seg.write(get_field_buffer(var_nz,'i',int_size))
        cfd_seg.write(get_field_buffer(var_stag_x,'d',dbl_size))
        cfd_seg.write(get_field_buffer(var_stag_y,'d',dbl_size))
        cfd_seg.write(get_field_buffer(var_stag_z,'d',dbl_size))       
        cfd_seg.write(get_field_buffer(var_min,'d',dbl_size))
        cfd_seg.write(get_field_buffer(var_max,'d',dbl_size))        
        cfd_seg.write(get_string_field_buffer(var_mesh_name,max_str_len))
        cfd_seg.write(get_string_field_buffer(var_mesh_class,max_str_len)) 
    
        cfd_seg.write(get_array_field_buffer(var_data,'%sd'))
        #################################################################################################################################################
    
    
    cfd.close()
    cfd_seg.close()
    
    print ' '
    print 'Size of ', cfd_path, ' is ', os.path.getsize(cfd_path)/(1024*1024.0), ' MB.'
    print 'Size of ', cfd_seg_path, ' is ', os.path.getsize(cfd_seg_path)/(1024*1024.0), ' MB.'
    print ' '
  


