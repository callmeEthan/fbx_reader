#region fbx.gml
function fbx_init()
{
	#macro FBX_NAME 0
	#macro FBX_TYPE 1
	#macro FBX_DATA 2
	#macro FBX_TREE 3

	globalvar FBX_TEMP_LISTS, FBX_TEMP_LISTS, FBX_SPRITE_POOL, FBX_KTIME, FBX_VERTEX_FORMAT;
	FBX_TEMP_LISTS  = ds_list_create();
	FBX_SPRITE_POOL = ds_map_create();
	FBX_KTIME       = 1/46186158000;

	vertex_format_begin();
	vertex_format_add_position_3d();
	vertex_format_add_normal();
	vertex_format_add_texcoord();
	vertex_format_add_color();          // tangent
	vertex_format_add_color();          // blend bones
	vertex_format_add_color();          // blend weights
	FBX_VERTEX_FORMAT = vertex_format_end();
	
	
	enum FBX_Model {uuid,name,type,skin,morphers,vertices,normals,uv,bone_weights,
	    index_list,face_list,vertex_buffer,limb,material,tangents,buffer,sizeof}

	enum FBX_Model_Data {uuid,name,type,loc,geom_loc,rot,rot_ofs,rot_piv,
	    pre_rot,pst_rot,rot_ord,rot_alt_mat,geom_rot,sca,sca_ofs,sca_piv,
	    geom_sca,rot_act,parent,bind_mat,bind_imt,deformer_index, sizeof
	    }
	enum FBX_RotOrd {XYZ,XZY,YZX,YXZ,ZXY,ZYX}
	enum FBX_Pose_Data {uuid,matrix,sizeof}
	enum FBX_Deformer_Data {uuid,name,type,indexes,weights,transform,
	    transform_link,transform_associate_model,bone_index,sizeof
	    }
	enum FBX_AnimCurveNode_Data {uuid,name,type,limb,transform,sizeof}
	enum FBX_AnimCurve_Data {uuid,name,keytime,keyvalue,sizeof}
	enum FBX_Skin_Data {uuid,name,type,deformers,sizeof}
	enum FBX_ALayer_Data {uuid,name,nodes,length,sizeof}
	enum FBX_AStack_Data {uuid,name,sizeof}
	enum FBX_Connection_Data {index,data,link,res_index,res_data,sizeof}

	enum FBX_Pool {geom, limbs, pose, deformer, anode, acurve, alayer, astack,
	    skin, material, attrib, ties, bone_mtx, mesh_mtx, anim_index,
	    anim_length, anim_time, anim_count, sizeof}
    
	enum FBX_Anim {loc, rot, sca, sizeof }
	enum FBX_Node {x,y,z,sizeof}
	enum FBX_Timeline {time,value,sizeof}
	enum FBX_PBR_Mat    {albedo,normal,emissive,arm,sizeof}
	enum FBX_BuffTypes {bind, skin, geom, anim, sizeof}
}
function fbx_clear()
{
	ds_list_destroy(FBX_TEMP_LISTS);
	var spr, key;
	key = ds_map_find_first(FBX_SPRITE_POOL);
	while ds_map_exists(FBX_SPRITE_POOL,key)
	{
	    spr = FBX_SPRITE_POOL[? key];
	    if sprite_exists(spr) then sprite_delete(spr);
	    key = ds_map_find_next(FBX_SPRITE_POOL,key);
	}
	vertex_format_delete(FBX_VERTEX_FORMAT);
}
function fbx_timeline_value_get(timeline, time)
{
	if array_length_1d(timeline) == 0 then return 0;
	var tims = timeline[0];
	var vals = timeline[1];
	if array_length_1d(tims) == 0 return 0;
	if array_length_1d(vals) == 0 return 0;
	var t1, t2;
	t1 = 0;
	t2 = 0;
	var len = array_length_1d(tims);
	var size = len - 1;
	var i = -1;
	var j = 0;
	while ++i < size
	{
	    t1 = tims[i];
	    j = (i + 1) % len;
	    t2 = tims[j];
	    if time >= t1 and time <= t2 then break;
	}
	if time > t2 then return vals[size];

	var part = ( time - t1 ) / (t2 - t1);
	var v1 = vals[i];
	var v2 = vals[j];
	return part * (v2 - v1) + v1;
}
function fbx_timeline_value_get_scale(timeline, time)
{
	if array_length_1d(timeline) == 0 then return 1;
	var tims = timeline[0];
	var vals = timeline[1];
	if array_length_1d(tims) == 0 return 0;
	if array_length_1d(vals) == 0 return 0;
	var t1, t2;
	t1 = 0;
	t2 = 0;
	var len = array_length_1d(tims);
	var size = len - 1;
	var i = -1;
	var j = 0;
	while ++i < size
	{
	    t1 = tims[i];
	    j = (i + 1) % len;
	    t2 = tims[j];
	    if time >= t1 and time <= t2 then break;
	}
	if time > t2 then return vals[size];

	var part = ( time - t1 ) / (t2 - t1);
	var v1 = vals[i];
	var v2 = vals[j];
	return part * (v2 - v1) + v1;
}
function fbx_timeline_value_360(timeline, time)
{
	if array_length_1d(timeline) == 0 then return 0;
	var tims = timeline[0];
	var vals = timeline[1];
	if array_length_1d(tims) == 0 return 0;
	if array_length_1d(vals) == 0 return 0;
	var t1, t2;
	t1 = 0;
	t2 = 0;
	var len = array_length_1d(tims);
	var size = len - 1;
	var i = -1;
	var j = 0;
	while ++i < size
	{
	    t1 = tims[i];
	    j = (i + 1) % len;
	    t2 = tims[j];
	    if time >= t1 and time <= t2 then break;
	}
	if time > t2 then return vals[size];

	var part = ( time - t1 ) / (t2 - t1);
	var v1 = vals[i];
	var v2 = vals[j];
	return part * angle_difference(v2,v1) + v1;
}
function fbx_timeline_length(timeline)
{
	if array_length_1d(timeline) == 0 then return 0;
	var time = timeline[0];
	var size = array_length_1d(time) - 1;
	var len = size < 0 ? 0: time[size];
	return len;
}
function fbx_anode_length(anode)
{
	var lx,ly,lz;
	if array_length_1d(anode) == 0 then return 0;
	lx = fbx_timeline_length(anode[0]);
	ly = fbx_timeline_length(anode[1]);
	lz = fbx_timeline_length(anode[2]);
	return max(lx,ly,lz);
}
function fbx_alayer_length()
{
	//!#import Blank
	var layr/*:FBX_ALayer_Data*/ = argument0;
	;
	var maxtime, nodes, i, len, lt,lr,ls,time,node;
	maxtime = 0;
	nodes = layr[FBX_ALayer_Data.nodes];
	i = -1;
	len = array_length_1d(nodes);
	while ++i < len
	{
	    node = nodes[i];
	    if not is_array(node) then continue;
	    lt = fbx_anode_length(node[0]);
	    lr = fbx_anode_length(node[1]);
	    ls = fbx_anode_length(node[2]);
	    time = max(lt,lr,ls);
	    maxtime = max(maxtime, time);
	}
	return maxtime;
}
function fbx_limb_get_mat()
{
	//!#import Blank
	var limb/*:FBX_Model_Data*/ = argument0;
	;
	var loc = limb[FBX_Model_Data.loc];
	var rot = limb[FBX_Model_Data.rot];
	var sca = limb[FBX_Model_Data.sca];

	return fbx_build_matrix(loc[0],loc[1],loc[2],-rot[0],-rot[1],-rot[2],sca[0],sca[1],sca[2] );
}
function fbx_limb_geom_mat()
{
	//!#import Blank
	var limb/*:FBX_Model_Data*/ = argument0;
	;
	var loc = limb[FBX_Model_Data.geom_loc];
	var rot = limb[FBX_Model_Data.geom_rot];
	var sca = limb[FBX_Model_Data.geom_sca];

	return fbx_build_matrix(loc[0],loc[1],loc[2],-rot[0],-rot[1],-rot[2],sca[0],sca[1],sca[2] );
}
function fbx_limb_pre_rot()
{
	//!#import Blank
	var limb/*:FBX_Model_Data*/ = argument0;
	;
	var rot = limb[FBX_Model_Data.pre_rot];

	return fbx_build_matrix(0,0,0,-rot[0],-rot[1],-rot[2],1,1,1 );
}
function fbx_build_bind_pose()
{
	/// @arg poses
	/// @arg limbs
	//!#import Blank
	var poses = argument0, limbs = argument1;
	;
	var i, j, len, limb/*:FBX_Model_Data*/, poses_len, pose/*:FBX_Pose_Data*/,
	    lid, pid;
	var len = array_length_1d(limbs);
	var i = -1;
	while ++i < len
	{
	    limb = limbs[i];
	    lid = limb[FBX_Model_Data.uuid];
	    poses_len = array_length_1d(poses);
	    j = -1;
	    while ++j < poses_len
	    {
	        pose = poses[j];
	        pid = pose[FBX_Pose_Data.uuid];
	        if lid != pid then continue;
	        var mat = pose[FBX_Pose_Data.matrix];
	        limb[@FBX_Model_Data.bind_mat] = mat;
	        limb[@FBX_Model_Data.bind_imt] = fbx_matrix_inverse( mat );
	        break;
	    }
	}
}
function fbx_build_bone_weight_list()
{
	/// @arg models
	// for each geometry skin make list of bone index
	// for each deformer add its index and weight to bone index list
	//!#import Blank
	var models = argument0;
	var length, i, m/*:FBX_Model*/, skin/*:FBX_Skin_Data*/,defs, defdat/*:FBX_Deformer_Data*/,
	    bidx, defs_len, j, indexes, weights, ind_size, k, ind, array, bw, 
	    morphers, bone_index;
	length = array_length_1d(models);
	i = -1;
	while ++i < length 
	{
	    m = models[i];
	    skin = m[FBX_Model.skin];
	    if skin == noone then continue;
	    defs = skin[FBX_Skin_Data.deformers];
	    morphers = m[FBX_Model.morphers];
	    bidx = fbx_ds_list_create();
	    m[@FBX_Model.bone_weights]  = bidx;
	    defs_len = array_length_1d(defs);
	    j = -1;
	    while ++j < defs_len
	    {
	        defdat = defs[j];
	        indexes = defdat[FBX_Deformer_Data.indexes];
	        weights = defdat[FBX_Deformer_Data.weights];
	        bone_index = defdat[FBX_Deformer_Data.bone_index];
	        morphers[@ bone_index ] = j;
	        if is_array(indexes) then continue;
	        ind_size = ds_list_size(indexes);
	        k = -1;
	        while ++k < ind_size
	        {
	            ind = indexes[|k];
	            array = bidx[| ind];
	            if not is_array(array)
	            {
	                array = [];
	                bidx[| ind] = array;
	            }
	            bw = [ j,weights[| k] ];
	            fbx_array_push(array,bw);
	        }
	    }
	}
}
function fbx_build_tangents_list(lpos, luvs, lnrm)
{
	var ltng, psize, i, vid1, vid2, vid3, tid1, tid2, tid3, v1, v2, v3,
	    uv1, uv2, uv3, n1, n2, n3, t1, t2, t3, tng;
	ltng = fbx_ds_list_create();
	psize = ds_list_size(lpos) / 9;
	i = -1;
	while ++i < psize
	{
	    vid1 = i * 9;
	    vid2 = vid1+3;
	    vid3 = vid2+3;
	    tid1 = i * 6;
	    tid2 = tid1 + 2;
	    tid3 = tid2 + 2;
	    v1 = [ lpos[| vid1], lpos[| vid1+1], lpos[| vid1+2] ];
	    v2 = [ lpos[| vid2], lpos[| vid2+1], lpos[| vid2+2] ];
	    v3 = [ lpos[| vid3], lpos[| vid3+1], lpos[| vid3+2] ];
	    uv1 = [ luvs[| tid1], luvs[| tid1 + 1] ];
	    uv2 = [ luvs[| tid2], luvs[| tid2 + 1] ];
	    uv3 = [ luvs[| tid3], luvs[| tid3 + 1] ];
	    tng = fbx_build_tangent(v1,v2,v3,uv1,uv2,uv3);
	    n1 = [ lnrm[| vid1], lnrm[| vid1+1], lnrm[| vid1+2] ];
	    n2 = [ lnrm[| vid2], lnrm[| vid2+1], lnrm[| vid2+2] ];
	    n3 = [ lnrm[| vid3], lnrm[| vid3+1], lnrm[| vid3+2] ];
	    t1 = fbx_tangent_nrm(tng, n1 );
	    t2 = fbx_tangent_nrm(tng, n2 );
	    t3 = fbx_tangent_nrm(tng, n3 );
	    ds_list_add( ltng, t1[0], t1[1], t1[2] );
	    ds_list_add( ltng, t2[0], t2[1], t2[2] );
	    ds_list_add( ltng, t3[0], t3[1], t3[2] );
	}
	return ltng;
}
function fbx_build_vertex_buffers()
{
	/// @arg models
	/// @arg limbs
	//!#import Blank
	var models = argument0, limbs = argument1;
	var mat, length, i, m/*:FBX_Model*/, b, verts, norms, uvs, bweight, indexes, skin,
	    size, j, vid, tid, nid, vv, nn, pid, bwt, bsz, bones, weigh, k, bw, bc, wc,
	    limb/*:FBX_Model_Data*/, tangs;
	length = array_length_1d(models);
	i = -1;
	while ++i < length 
	{
	    m = models[i];
	    // here should be triangulation
	    verts   = m[FBX_Model.vertices];
	    norms   = m[FBX_Model.normals];
	    tangs   = m[FBX_Model.tangents];
	    uvs     = m[FBX_Model.uv];
	    bweight = m[FBX_Model.bone_weights];
	    indexes = m[FBX_Model.index_list];
	    skin    = m[FBX_Model.skin];
	    limb    = limbs[ m[FBX_Model.limb] ];
	    mat     = limb[FBX_Model_Data.bind_mat];
	    size = ds_list_size(verts) / 3;
	    if size == 0
	    {
	        continue;
	    }
    
	    tangs = fbx_build_tangents_list(verts, uvs, norms );
    
	    b = vertex_create_buffer();
	    vertex_begin(b, FBX_VERTEX_FORMAT);
	    j = -1;
	    if skin == noone
	    {
	        while ++j < size
	        {
	            vid = j * 3;
	            tid = j * 2;
	            nid = vid;
	            vertex_position_3d(b,verts[|vid],verts[|vid+1],verts[|vid+2]);
	            vertex_normal(b,norms[|nid],norms[|nid+1],norms[|nid+2]);
	            vertex_texcoord(b,uvs[|tid],1-uvs[|tid+1]);
	            var tang = [tangs[|nid],tangs[|nid+1],tangs[|nid+2]];
	            tang = fbx_convert_vector_to_rgb(tang);
	            vertex_color(b,tang,1);
	            vertex_color(b,0,0);
	            vertex_color(b,c_white,1);
	        }
	    }
	    else 
	    {
	        while ++j < size
	        {
	            vid = j * 3;
	            tid = j * 2;
	            nid = vid;
	            vv = [ verts[|vid],verts[|vid+1],verts[|vid+2] ];
	            nn = [ norms[|nid],norms[|nid+1],norms[|nid+2] ];
	            vv = matrix_transform_vertex(mat,vv[0],vv[1],vv[2]);
	            nn = matrix_transform_vertex(mat,nn[0],nn[1],nn[2]);
	            vertex_position_3d(b,vv[0],vv[1],vv[2]);
	            vertex_normal(b,nn[0],nn[1],nn[2]);
	            vertex_texcoord(b,uvs[|tid],1-uvs[|tid+1]);
	            var tang = [tangs[|nid],tangs[|nid+1],tangs[|nid+2]];
	            tang = fbx_convert_vector_to_rgb(tang);
	            vertex_color(b,tang,1);
            
	            pid = indexes[| j];
	            bwt = bweight[| pid];
	            bsz = array_length_1d(bwt);
	            bones = [];
	            weigh = [];
	            k = -1;
	            while ++k < 4
	            {
	                if k < bsz
	                {
	                    bw = bwt[k];
	                    bones[k] = bw[0];
	                    weigh[k] = bw[1] * 255;
	                }
	                else
	                {
	                    bones[k] = 0;
	                    weigh[k] = 0;
	                }
	            }
            
	            bc = fbx_color_array_to_rgb( bones );
	            wc = fbx_color_array_to_rgb( weigh );
            
	            vertex_color(b,bc,bones[3]);
	            vertex_color(b,wc,weigh[3]);
	        }
	    }
	    vertex_end(b);
	    m[@FBX_Model.buffer] = buffer_create_from_vertex_buffer(b, buffer_fixed, 1 );
	    vertex_freeze(b);
	    m[@FBX_Model.vertex_buffer] = b;
	}
}
function fbx_ds_list_create()
{
	var list = ds_list_create();
	fbx_ds_list_add_list(FBX_TEMP_LISTS, list);
	return list;
}
function fbx_draw()
{
	//!#import blank.* as Blank
	var d/*:FBX_Pool*/ = argument0;
	var transmat = d[FBX_Pool.mesh_mtx];
	var models = d[FBX_Pool.geom];
	var bone_mtx = d[FBX_Pool.bone_mtx];
	var sample = [];
	var world_mat = matrix_get(matrix_world);

	var len = array_length_1d(models);
	var i = -1;
	while ++i < len
	{
	    var m/*:FBX_Model*/ = models[i];
	    var b = m[FBX_Model.vertex_buffer];
	    var s/*:FBX_Skin_Data*/ = m[FBX_Model.skin];
	    sample = [];
	    if s == noone
	    {
	        sample = matrix_build_identity();
	    }
	    else
	    {
	        var morphers = s[FBX_Skin_Data.deformers];
	        var mlen = array_length_1d(morphers);
	        sample = [];
	        var j = -1;
	        while ++j < mlen
	        {
	            var morpher/*:FBX_Deformer_Data*/ = morphers[j];
	            var bindex = morpher[FBX_Deformer_Data.bone_index];
	            fbx_array_concat(sample, transmat[bindex] );
	        }
	    }
    
	    var limb_index = m[FBX_Model.limb];
	    matrix_set(matrix_world, matrix_multiply( bone_mtx[limb_index], world_mat ) );
    
	    var u_bones = shader_get_uniform(shader_current(), "u_bones");
	    shader_set_uniform_f_array( u_bones, sample);
	    var albedo = m[FBX_Model.material];
	    if is_array(albedo)
	    {
	        albedo = fbx_material_set(albedo);
	    }
	    if b != noone 
	    {
	        vertex_submit(b,pr_trianglelist,albedo);
	    }
	}
}
function fbx_frame_matrix(anim, time)
{
	var atx,aty,atz, arx,ary,arz, asx,asy,asz, at, ar, as;
	at  = anim[0];
	if array_length_1d(at) == 0 or !is_array(at)
	{
	    atx = 0;aty = 0;atz = 0;
	}
	else
	{
	    atx = fbx_timeline_value_get(at[0],time);
	    aty = fbx_timeline_value_get(at[1],time);
	    atz = fbx_timeline_value_get(at[2],time);
	}
	ar  = anim[1];
	if array_length_1d(ar) == 0 or !is_array(ar)
	{
	    arx = 0;ary = 0;arz = 0;
	}
	else
	{
	    arx = -fbx_timeline_value_360(ar[0],time);
	    ary = -fbx_timeline_value_360(ar[1],time);
	    arz = -fbx_timeline_value_360(ar[2],time);
	}
	as  = anim[2];
	if array_length_1d(as) == 0 then {asx = 1;asy = 1;asz = 1;}
	else
	{
	    asx = fbx_timeline_value_get_scale(as[0],time);
	    asy = fbx_timeline_value_get_scale(as[1],time);
	    asz = fbx_timeline_value_get_scale(as[2],time);
	}

	return fbx_build_matrix(atx,aty,atz,arx,ary,arz,asx,asy,asz);
}

function fbx_animate()
{
	//!#import Blank
	var fbx/*:FBX_Pool*/ = argument0, tick = argument1;
	;
	var limbs   = fbx[FBX_Pool.limbs];
	var alayer = fbx[FBX_Pool.alayer];
	if fbx[FBX_Pool.anim_count] == 0 return false;
	var anim_index = fbx[FBX_Pool.anim_index];
	var lay/*:FBX_ALayer_Data*/ = alayer[anim_index];
	var limb_anim = lay[FBX_ALayer_Data.nodes];
	var bone_mtx = fbx[FBX_Pool.bone_mtx];
	var mesh_mtx = fbx[FBX_Pool.mesh_mtx];
	var size = array_length_1d(limb_anim);
	var time = fbx[FBX_Pool.anim_time] + tick;
	var anim_len = fbx[FBX_Pool.anim_length];
	if time >= anim_len then time -= anim_len;
	fbx[@FBX_Pool.anim_time] = time;

	var i = -1;
	while ++i < size
	{
	    var anim = limb_anim[i];
	    if not is_array(anim) then continue;
	    var amt = fbx_frame_matrix( anim, time );

	    var limb/*:FBX_Model_Data*/ = limbs[i];
	    var parent = limb[FBX_Model_Data.parent];
	    if parent != noone
	    {
	        amt = matrix_multiply(amt,bone_mtx[parent]);
	    }
	    bone_mtx[@ i] = amt;
    
	    mesh_mtx[@ i] = matrix_multiply(limb[FBX_Model_Data.bind_imt],amt);
	}
}
function fbx_read(fbx)
{
	/// @arg fbx
	//!#import Blank
	log("fbx: Collecting data...")
	var d/*:FBX_Pool*/ = fbx_collect_data( fbx );

	log("fbx: Connecting data...")
	// here we use script names instead of sctipt indexes
	// because if script defined in extension then it is possible to
	// get correct index for it only using asset_get_index
	fbx_connect( "fbx_connect_curve_to_node", d[FBX_Pool.anode],d[FBX_Pool.acurve],d[FBX_Pool.ties]);
	fbx_connect( "fbx_connect_limbs", d[FBX_Pool.limbs],d[FBX_Pool.limbs],d[FBX_Pool.ties]);
	fbx_connect( "fbx_connect_bones_to_deformers", d[FBX_Pool.deformer],d[FBX_Pool.limbs], d[FBX_Pool.ties] );
	fbx_connect( "fbx_connect_clusters_to_skins", d[FBX_Pool.skin], d[FBX_Pool.deformer], d[FBX_Pool.ties]);
	fbx_connect( "fbx_connect_anode_to_limbs", d[FBX_Pool.limbs], d[FBX_Pool.anode], d[FBX_Pool.ties]);
	fbx_connect( "fbx_connect_anode_to_alayer", d[FBX_Pool.alayer], d[FBX_Pool.anode], d[FBX_Pool.ties]);
	fbx_connect( "fbx_connect_skin_to_geometry", d[FBX_Pool.geom], d[FBX_Pool.skin], d[FBX_Pool.ties]);
	fbx_connect( "fbx_connect_geom_to_limb",d[FBX_Pool.limbs],d[FBX_Pool.geom],d[FBX_Pool.ties]);

	log("fbx: Binding pose...")
	fbx_build_bind_pose(d[FBX_Pool.pose],d[FBX_Pool.limbs]);
	log("fbx: Bone weight list...")
	fbx_build_bone_weight_list(d[FBX_Pool.geom]);
	log("fbx: Building vertex buffers...")
	fbx_build_vertex_buffers(d[FBX_Pool.geom],d[FBX_Pool.limbs]);

	/// BUILD BONES
	log("fbx: Building bones...")
	var bone_matrix = [];
	var mesh_matrix = [];
	var model_data = d[FBX_Pool.limbs];
	var size, i, md/*:FBX_Model_Data*/, par_id, mtx, inst;
	var size = array_length_1d(model_data);
	var i = -1;
	while ++i < size
	{
	    var md/*:FBX_Model_Data*/ = model_data[i];
	    mtx = md[FBX_Model_Data.bind_mat];
	    // var lmb_mat = fbx_limb_get_mat( md );
	    // var pre_rot = fbx_limb_pre_rot(md);
	    // var geo_mat = fbx_limb_geom_mat( md );
	    // pre_rot = matrix_multiply( pre_rot, geo_mat );
	    // if i == 0 then mtx = matrix_multiply( lmb_mat, pre_rot );
	    // else mtx = matrix_multiply( pre_rot, lmb_mat );
	    // var inst = instance_create_layer(0,0,layer,obj_bone);
	    // inst.mtx = mtx;
	    // inst.index = i;
	    bone_matrix[i] = matrix_build_identity();
	    mesh_matrix[i] = matrix_build_identity();
	}
	log("fbx: Tracing number of bones...")
	fbx_trace( "NUMBER OF BONES: ", size );

	d[@FBX_Pool.bone_mtx] = bone_matrix;
	d[@FBX_Pool.mesh_mtx] = mesh_matrix;

	var alayer = d[FBX_Pool.alayer];
	var alength = array_length_1d( alayer );
	d[@FBX_Pool.anim_count] = alength;
	d[@FBX_Pool.anim_time] = 0;

	var i = -1;
	while ++i < alength
	{
	    var anim_index = d[FBX_Pool.anim_index];
	    var lay/*:FBX_ALayer_Data*/ = alayer[anim_index];
	    lay[@FBX_ALayer_Data.length] = fbx_alayer_length(lay);
	}

	log("fbx: Finalizing...")
	ds_list_clear(FBX_TEMP_LISTS);
	ds_list_destroy(fbx);

	return d;
}
function fbx_animation_set()
{
	//!#import Blank
	var fbx/*:FBX_Pool*/ = argument0, index = argument1;
	var anim_count, layers, lay;
	anim_count = fbx[FBX_Pool.anim_count];
	if anim_count == 0 then return 0;
	index = index % anim_count;
	var layers = fbx[FBX_Pool.alayer];
	var lay/*:FBX_ALayer_Data*/ = layers[index];
	fbx[@FBX_Pool.anim_index] = index;
	fbx[@FBX_Pool.anim_time] = 0;
	fbx[@FBX_Pool.anim_length] = lay[FBX_ALayer_Data.length];
}
function fbx_load(filename)
{
	var fbx = fbx_parse(filename);

	return fbx_read(fbx);
}
function fbx_save_json(fbx)
{
	var map = ds_map_create();
	ds_map_add_list(map,"FBX",fbx);
	var jencoded = json_encode(map);
	var file = file_text_open_write("fbx.json");
	file_text_write_string(file,jencoded);
	file_text_close(file);
}
#endregion

#region fbx_parse.gml
function fbx_parse(filename)
{
	var b = buffer_load(filename);
	var head_magic = fbx_array_from_string("Kaydara FBX Binary");
	//fbx_array_concat(head_magic,[$20,$20,$00,$1A,$00]);
	fbx_array_concat(head_magic,[32,32,0,26,0]);

	var magic = fbx_buffer_read_len(b,array_length_1d(head_magic));
	if not array_equals(head_magic,magic) 
	{
	    show_error( "FBX Invalid Header", true );
	    return noone;
	}

	var version             = buffer_read(b,buffer_u32);
	show_debug_message( "fbx_version_is " + string(version) );
	if version < 7100
	{
	    show_error( "FBX version older then 7100 is not supported", true );
	    return noone;
	}

	var block_len   = 25;
	var elem_type   = buffer_u64;
	if version < 7500 
	{
	    block_len   = 13;
	    elem_type   = buffer_u32;
	}

	var elements = ds_list_create();
	while true 
	{
	    var elem = fbx_parse_read_elem( b, elem_type, block_len );
	    if elem == noone then break;
	    fbx_ds_list_add_list(elements,elem);
	}

	var root = ds_list_create();
	ds_list_add(root,"root");
	ds_list_add(root,"");
	fbx_ds_list_add_list(root,ds_list_create());
	fbx_ds_list_add_list(root,elements);

	buffer_delete(b);
	return root;
}

function fbx_buffer_read_string_ubyte(buffer)
{
	var size, array;
	size = buffer_read(buffer,buffer_u8);
	array= fbx_buffer_read_len(buffer,size);
	return fbx_string_from_array( array );
}

function fbx_parse_read_array(buffer, data_type)
{
	var length, encoding, comp_len, data, stride, decompressed, list, i;
	length = buffer_read(buffer,buffer_u32);
	encoding = buffer_read(buffer,buffer_u32);
	comp_len = buffer_read(buffer,buffer_u32);
	if comp_len == 0 
	{
	    return ds_list_create();
	}
	data = fbx_buffer_read_data(buffer,comp_len);
	stride = buffer_sizeof(data_type);
	if encoding == 1 {
	    decompressed = buffer_decompress(data);
	    buffer_delete(data);
	    data = decompressed;}

	fbx_assert(length * stride == buffer_get_size(data), "Wrong size of decompressed buffer" )

	list = ds_list_create();
	i = -1;
	while ++i < length
	    ds_list_add(list,buffer_read(data,data_type));
	buffer_delete(data);
	return list;
}

function fbx_parse_read_type(buffer, type)
{
	var size, array, data;
	switch type 
	{
	    case "Y": return buffer_read(buffer,buffer_s16);break;
	    case "C": return buffer_read(buffer,buffer_u8);break;
	    case "I": return buffer_read(buffer,buffer_s32);break;
	    case "F": return buffer_read(buffer,buffer_f32);break;
	    case "D": return buffer_read(buffer,buffer_f64);break;
	    case "L": return buffer_read(buffer,buffer_u64);break; /// here should be s64, but it seems u64 actually s64
	    case "R": 
	        size = buffer_read(buffer,buffer_u32);
	        array = fbx_buffer_read_len(buffer,size);
	        return fbx_ds_list_from_array( array );
	    break;
	    case "S":
	        size = buffer_read(buffer,buffer_u32);
	        data = fbx_buffer_read_len(buffer,size);
	        return fbx_string_from_array(data);
	    break;
	    case "f": return fbx_parse_read_array( buffer, buffer_f32);break;
	    case "i": return fbx_parse_read_array( buffer, buffer_s32);break;
	    case "d": return fbx_parse_read_array( buffer, buffer_f64);break;
	    case "l": return fbx_parse_read_array( buffer, buffer_u64);break; /// here should be s64, but it seems u64 actually s64
	    case "buffer": return fbx_parse_read_array( buffer, buffer_u8);break;
	    case "c": return fbx_parse_read_array( buffer, buffer_u8);break;
	}
}

function fbx_parse_read_elem(b, element_type, block_len)
{
	var end_offset      = buffer_read(b,element_type);
	if end_offset  == 0 return noone;

	var prop_count  = buffer_read(b,element_type);
	var prop_length = buffer_read(b,element_type);
	var prop_name   = fbx_buffer_read_string_ubyte(b);
	var props_type  = "";
	var props_data  = ds_list_create();
	var elem_tree   = ds_list_create();

	var i = -1;

	while ++i < prop_count
	{
	    var data_type = chr(buffer_read(b,buffer_u8));
	    var data = fbx_parse_read_type(b,data_type);
	    props_data[| i] = data;
	    var list_types = "Rfidlbc";
	    var is_list = string_pos(data_type,list_types) > 0;
	    if is_list then ds_list_mark_as_list(props_data,i);
	    props_type += data_type;
	}

	var tell = buffer_tell(b);
	if buffer_tell(b) < end_offset 
	{
	    while buffer_tell(b) < (end_offset - block_len )
	    {
	        var elem = fbx_parse_read_elem(b,element_type,block_len);
	        fbx_ds_list_add_list(elem_tree, elem);
	    }
    
	    var block_sentinel = fbx_buffer_read_len( b, block_len );
	    if not array_equals(block_sentinel,array_create(block_len,0) )
	    {
	        show_error( "IOError " + string( buffer_tell(b) ) + 
	            " failed to read nested block sentinel, expected all bytes to be 0", 
	            false);
	    }
	}

	if buffer_tell(b) != end_offset
	{
	    show_error( "FBX Scope length not reached, something is wrong. ", false );
	}

	var args = ds_list_create();
	ds_list_add(args,prop_name);
	ds_list_add(args,props_type);
	fbx_ds_list_add_list(args,props_data);
	fbx_ds_list_add_list(args,elem_tree);

	return args;
}
#endregion

#region fbx_common.gml
function fbx_array_push(array, value)
{
var len = array_length( array );
array[@len] = value;
}

function fbx_array_from_string(str)
{
	var i, array;
	array = [];
	i = -1;
	while ++i < string_length(str)
	    fbx_array_push(array,string_byte_at(str,i+1));
	return array;
}

function fbx_array_concat(a1, a2)
{
	var i = -1;
	var len = array_length(a2);
	while ++i < len
	{
	    fbx_array_push(a1,a2[i]);
	}
}

function fbx_array_from_list(list)
{
	var a, size, i;
	a = [];
	size = ds_list_size(list);
	i = -1;
	while ++i < size
	{
	    fbx_array_push(a,list[| i]);
	}
	return a;
}

function fbx_buffer_read_len(buff, len)
{
	var array, byte;
	array   = [];
	repeat len
	{
	    var byte = buffer_read(buff, buffer_u8);
	    fbx_array_push(array,byte);
	}
	return array;
}

function fbx_buffer_read_data(buffer, size)
{
	var out, offset;
	out = buffer_create(size,buffer_fixed,1);
	offset = buffer_tell(buffer);
	buffer_copy(buffer,offset,size,out,0);
	buffer_seek(buffer,buffer_seek_start,offset + size);
	return out;
}

function fbx_ds_list_add_list(list, value)
{
	ds_list_add(list,value);
	ds_list_mark_as_list(list,ds_list_size(list)-1);
}

function fbx_ds_list_add_map(list, map)
{
	ds_list_add(list,map);
	ds_list_mark_as_map(list,ds_list_size(list)-1);
}

function fbx_ds_list_from_array(array)
{
	var len, list, i;
	len = array_length(array);
	list = ds_list_create();
	i = -1;
	while ++i < len
	    list[|i] = array[i];
	return list;
}

function fbx_string_from_array(array)
{
	var i, str, len;
	i   = -1;
	str = "";
	len = array_length( array );
	while ++i < len
	{
	    str += chr(array[i]);
	}
	return str;
}

function fbx_string_to_array(str, delim)
{
	var pos, ind, a, len;
	pos = string_pos(delim,str);
	ind = 0;
	a = [];
	while(pos > 0) {
	    var s = string_copy(str,1,pos - 1);
	    a[ind] = s;
	    ind++;
	    len = string_length(str) - pos;
	    str = string_copy(str,pos + 1,len);
	    pos = string_pos(delim,str);
	}
	a[ind] = string_copy(str,1,string_length(str));
	return a;
}


function fbx_assert()
{
	if argument0 then exit;
	show_error("FALSE_ASSERT " + string( argument1 ), true );
}
function fbx_trace()
{
	var str = string(get_timer()) + " ";
	var i = -1;
	while(++i < argument_count) str += string(argument[i]) + " ";
	show_debug_message(str);
}

function fbx_convert_vector_to_rgb(v)
{
	var vx,vy,vz;
	vx = fbx_convert_value_to_rgb( v[0] );
	vy = fbx_convert_value_to_rgb( v[1] );
	vz = fbx_convert_value_to_rgb( v[2] );
	return fbx_color_array_to_rgb( [vx,vy,vz] );
}
function fbx_convert_value_to_rgb(value)
{
	return round( (value + 1) * 127.5 );
}

function fbx_color_array_to_rgb(clr)
{
	var rgb;
	rgb = clr[0]<<16;
	rgb += clr[1]<<8;
	rgb += clr[2];
	return rgb;
}

function fbx_build_tangent()
{
	/// @arg p1
	/// @arg p2
	/// @arg p3
	/// @arg t1
	/// @arg t2
	/// @arg t3
	var p1 = argument0, p2 = argument1, p3 = argument2, t1 = argument3, t2 = argument4, t3 = argument5;
	var dpx, dpy, du, dv, duv, tng;
	dpx = fbx_v3_sub( p2, p1 );
	dpy = fbx_v3_sub( p3, p1 );
	du = fbx_v2_sub( t2, t1 );
	dv = fbx_v2_sub( t3, t1 );
	duv = du[0] * dv[1] - dv[0] * du[1];
	duv = duv == 0?0:1 / duv;
	dpx = fbx_v3_mul( dpx, dv[1] );
	dpy = fbx_v3_mul( dpy, du[1] );
	tng = fbx_v3_sub( dpx, dpy );
	return fbx_v3_mul( tng, duv );
}

function fbx_tangent_nrm(tng, nrm)
{
	var dp;
	nrm = fbx_v3_nrm( nrm );
	dp = fbx_v3_dot( nrm, tng );
	nrm = fbx_v3_mul( nrm, dp );
	tng = fbx_v3_sub( tng, nrm );
	tng = fbx_v3_nrm( tng );
	return tng;
}

function fbx_v2_sub()
{
	///fbx_v2_add(a:v2,b:v2)->v2 :subtract vector b from vector a
	return [argument0[0] - argument1[0], argument0[1] - argument1[1]];
}
function fbx_v3_mul()
{
	///fbx_v3_mul(a:v3,m:real)->v3 :multiply vector by number
	return [argument0[0] * argument1, argument0[1] * argument1, argument0[2] * argument1];
}
function fbx_v3_nrm()
{
	///fbx_v3_nrm(a:v3)->v3 :normalize vector
	var l = point_distance_3d(0,0,0,argument0[0],argument0[1],argument0[2]);
	l = l==0?0:1/l;
	return [argument0[0] * l, argument0[1] * l, argument0[2] * l];
}

function fbx_v3_dot()
{
	///fbx_v3_dot(a:v3,b:v3)->real :return dot product
	return dot_product_3d(argument0[0],argument0[1],argument0[2],argument1[0],argument1[1],argument1[2]);
}
function fbx_v3_sub()
{
	///fbx_v3_sub(a:v3,b:v3)->v3 :subtract vector b from vector a
	return [argument0[0] - argument1[0], argument0[1] - argument1[1], argument0[2] - argument1[2]];
}

function fbx_build_matrix(px, py, pz, rx, ry, rz, sx, sy, sz)
{
	var mrx = matrix_build(0,0,0,rx,0,0,1,1,1);
	var mry = matrix_build(0,0,0,0,ry,0,1,1,1);
	var mrz = matrix_build(0,0,0,0,0,rz,1,1,1);
	var mtr = matrix_build(px,py,pz,0,0,0,sx,sy,sz);
	var amt = matrix_multiply( matrix_multiply(mrx,mry),mrz );
	amt = matrix_multiply(amt,mtr);
	return amt;
}

function fbx_matrix_inverse(m)
{
	var a1 = m[0], a2 = m[1], a3 = m[2], a4 = m[3],
	    b1 = m[4], b2 = m[5], b3 = m[6], b4 = m[7],
	    c1 = m[8], c2 = m[9], c3 = m[10], c4 = m[11],
	    d1 = m[12], d2 = m[13], d3 = m[14], d4 = m[15],
	    i = -1, det, inv,
	    cd34=c3*d4-c4*d3, bd34=b3*d4-b4*d3, bc34=b3*c4-b4*c3,
	    ad34=a3*d4-a4*d3, ab34=a3*b4-a4*b3, ac34=a3*c4-a4*c3,
	    cd24=c2*d4-c4*d2, ad24=a2*d4-a4*d2, bc24=b2*c4-b4*c2,
	    bd24=b2*d4-b4*d2, ab24=a2*b4-a4*b2, ac24=a2*c4-a4*c2,
	    cd23=c2*d3-c3*d2, bd23=b2*d3-b3*d2, ad23=a2*d3-a3*d2,
	    ac23=a2*c3-a3*c2, bc23=b2*c3-b3*c2, ab23=a2*b3-a3*b2;
	inv = [
	     b2*cd34 - c2*bd34 + d2*bc34,
	    -a2*cd34 + c2*ad34 - d2*ac34,
	     a2*bd34 - b2*ad34 + d2*ab34,
	    -a2*bc34 + b2*ac34 - c2*ab34,
    
	    -b1*cd34 + c1*bd34 - d1*bc34,
	     a1*cd34 - c1*ad34 + d1*ac34,
	    -a1*bd34 + b1*ad34 - d1*ab34,
	     a1*bc34 - b1*ac34 + c1*ab34,
    
	     b1*cd24 - c1*bd24 + d1*bc24,
	    -a1*cd24 + c1*ad24 - d1*ac24,
	     a1*bd24 - b1*ad24 + d1*ab24,
	    -a1*bc24 + b1*ac24 - c1*ab24,
    
	    -b1*cd23 + c1*bd23 - d1*bc23,
	     a1*cd23 - c1*ad23 + d1*ac23,
	    -a1*bd23 + b1*ad23 - d1*ab23,
	     a1*bc23 - b1*ac23 + c1*ab23
	];

	det = a1*inv[0] + a2*inv[4] + a3*inv[8] + a4*inv[12];
	if det == 0 return m;
	det = 1 / det;
	while ++i < 16
	{
	    inv[i] *= det;
	}
	return inv;
}

function fbx_sprite_load(name)
{
	if ds_map_exists(FBX_SPRITE_POOL,name)
	{
	    return FBX_SPRITE_POOL[? name];
	}
	else
	{
	    var spr = sprite_add(name,0,false,true,0,0);
	    FBX_SPRITE_POOL[? name] = spr;
	    return spr;
	}
}

function fbx_material_set()
{
	/// @arg mat
	//!#import Blank
	var mat/*:FBX_PBR_Mat*/ = argument0;
	;
	fbx_shader_set_sampler( "u_normal", mat[FBX_PBR_Mat.normal] );
	fbx_shader_set_sampler( "u_emission", mat[FBX_PBR_Mat.emissive] );
	fbx_shader_set_sampler( "u_MetalicRoughnessSampler", mat[FBX_PBR_Mat.arm] );
	return mat[FBX_PBR_Mat.albedo];
}

function fbx_shader_set_sampler()
{
	/// @arg tex
/// @arg name
var name = argument0, tex = argument1;
var sampler = shader_get_sampler_index( shader_current(), name );
texture_set_stage( sampler, tex );
}

function fbx_material_assign()
{
	/// @arg fbx
	/// @arg name
	//! #import Blank
	var fbx/*:FBX_Pool*/ = argument0, name = argument1;
	;
	var geoms, mat, i, len, geom/*:FBX_Model*/;
	mat     = fbx_load_pbr_mat(name);
	geoms   = fbx[FBX_Pool.geom];
	len     = array_length(geoms);
	i = -1;
	while ++i < len
	{
	    geom = geoms[i];
	    geom[@FBX_Model.material] = mat;
	}
}

function fbx_material_assign_by_index()
{
	//! #import Blank
	var fbx/*:FBX_Pool*/ = argument0, index = argument1, name = argument2;
	;
	var geoms, mat, i, len, geom/*:FBX_Model*/;
	mat     = fbx_load_pbr_mat(name);
	geoms   = fbx[FBX_Pool.geom];
	len     = array_length(geoms);
	i = -1;
	while ++i < len
	{
	    if i == index
	    {
	        geom = geoms[i];
	        geom[@FBX_Model.material] = mat;
	        break;
	    }
	}
}

function fbx_load_pbr_mat(name)
{
	var spr_alb, spr_arm, spr_nrm, spr_emv, mat/*:FBX_PBR_Mat*/;
	var prefix = "assets/";
	spr_alb     = fbx_sprite_load( prefix + name + "_alb.jpg" );
	spr_arm     = fbx_sprite_load( prefix + name + "_arm.jpg" );
	spr_nrm     = fbx_sprite_load( prefix + name + "_nrm.jpg" );
	spr_emv     = fbx_sprite_load( prefix + name + "_emv.jpg" );
	mat = array_create(FBX_PBR_Mat.sizeof);
	mat[@FBX_PBR_Mat.albedo]      = sprite_get_texture( spr_alb, 0 );
	mat[@FBX_PBR_Mat.normal]      = sprite_get_texture( spr_nrm, 0 );
	mat[@FBX_PBR_Mat.emissive]    = sprite_get_texture( spr_emv, 0 );
	mat[@FBX_PBR_Mat.arm]         = sprite_get_texture( spr_arm, 0 );
	return mat;
}
#endregion

#region fbx_collect.gml
function fbx_collect_model_data()
{
	/// @arg node
	/// @arg array
	//!#import fbx_elem.* as FBX_Elem
	var node/*:FBX_Elem*/ = argument0, array = argument1;
	var dat/*:FBX_Model_Data*/ = array_create(FBX_Model_Data.sizeof);
	dat[@FBX_Model_Data.uuid] = fbx_elem_uuid(node);
	dat[@FBX_Model_Data.name] = fbx_elem_name(node);
	dat[@FBX_Model_Data.type] = fbx_elem_prop_value(node, 2);
	dat[@FBX_Model_Data.loc] = fbx_elem_prop_vector(node, "Lcl Translation",[0,0,0]);
	dat[@FBX_Model_Data.rot] = fbx_elem_prop_vector(node, "Lcl Rotation",[0,0,0]);
	dat[@FBX_Model_Data.sca] = fbx_elem_prop_vector(node, "Lcl Scaling",[1,1,1]);
	dat[@FBX_Model_Data.geom_loc] = fbx_elem_prop_vector(node, "GeometricTranslation",[0,0,0]);
	dat[@FBX_Model_Data.geom_rot] = fbx_elem_prop_vector(node, "GeometricRotation",[0,0,0]);
	dat[@FBX_Model_Data.geom_sca] = fbx_elem_prop_vector(node, "GeometricScaling",[1,1,1]);
	dat[@FBX_Model_Data.rot_ofs] = fbx_elem_prop_vector(node, "RotationOffset",[0,0,0]);
	dat[@FBX_Model_Data.rot_piv] = fbx_elem_prop_vector(node, "RotationPivot",[0,0,0]);
	dat[@FBX_Model_Data.sca_ofs] = fbx_elem_prop_vector(node, "ScalingOffset",[0,0,0]);
	dat[@FBX_Model_Data.sca_piv] = fbx_elem_prop_vector(node, "ScalingPivot",[0,0,0]);
	dat[@FBX_Model_Data.pre_rot] = fbx_elem_prop_vector(node, "PreRotation",[0,0,0]);
	dat[@FBX_Model_Data.pst_rot] = fbx_elem_prop_vector(node, "PostRotation",[0,0,0]);
	dat[@FBX_Model_Data.rot_act] = fbx_elem_prop_bool(node, "RotationActive",false);
	dat[@FBX_Model_Data.rot_ord] = fbx_elem_prop_enum(node, "RotationOrder",0);
	dat[@FBX_Model_Data.bind_mat] = matrix_build_identity();
	dat[@FBX_Model_Data.bind_imt] = matrix_build_identity();
	dat[@FBX_Model_Data.parent]  = noone;
	dat[@FBX_Model_Data.deformer_index] = noone;
	fbx_array_push(array, dat);
}

function fbx_collect_pose_data()
{
	/// @arg pose
	/// @arg array
	//!#import fbx_elem.* as FBX_Elem
	var pose/*:FBX_Elem*/ = argument0, array = argument1;
	var tree = fbx_elem_elems(pose);
	var length = ds_list_size(tree);
	var i = -1;
	while ++i < length 
	{
	    var node/*:FBX_Elem*/ = tree[|i];
	    if fbx_elem_id(node) != "PoseNode" then continue;
	    var nd/*:FBX_Elem*/ = fbx_elem_find_first(node, "Node");
	    var dat/*:FBX_Pose_Data*/ = array_create(FBX_Pose_Data.sizeof);
	    dat[@FBX_Pose_Data.uuid] = fbx_elem_prop_value(nd, 0);
	    dat[@FBX_Pose_Data.matrix] = fbx_array_from_list( fbx_elem_property(node, "Matrix",0,[]) );
	    fbx_array_push(array,dat);
	}
}

function fbx_collect_geometry_data()
{
	/// @arg node
	/// @arg array
	//!#import fbx_elem.* as FBX_Elem
	var node/*:FBX_Elem*/ = argument0, array = argument1;
	var vert_list, face_list, index_list, vertices, uv, normals, mdl/*:FBX_Model*/,
	    tangents;
	vert_list           = fbx_elem_property(node, "Vertices",0,-4);
	face_list           = fbx_elem_property(node, "PolygonVertexIndex",0,-4 );

	index_list      = fbx_ds_list_create();
	vertices        = fbx_create_vertices_list(vert_list,face_list, index_list);
	uv              = fbx_elem_geom_layer(node, "LayerElementUV", index_list);
	normals         = fbx_elem_geom_layer(node, "LayerElementNormal", index_list);
	tangents        = fbx_elem_geom_layer(node, "LayerElementTangent", index_list );
	// trace( node.name(), ds_list_size( vertices ), 
	//     ds_list_size( normals ), ds_list_size( uv ) );
	if uv == noone or normals == noone then return false;
	mdl = array_create(FBX_Model.sizeof);
	mdl[@FBX_Model.uuid]        = fbx_elem_uuid(node);
	mdl[@FBX_Model.name]        = fbx_elem_name(node);
	mdl[@FBX_Model.type]        = fbx_elem_prop_value(node, 2);
	mdl[@FBX_Model.skin]        = noone;
	mdl[@FBX_Model.vertices]    = vertices;
	mdl[@FBX_Model.normals]     = normals;
	mdl[@FBX_Model.tangents]    = tangents;
	mdl[@FBX_Model.uv]          = uv;
	mdl[@FBX_Model.index_list]  = index_list;
	mdl[@FBX_Model.face_list]   = face_list;
	mdl[@FBX_Model.morphers]    = [];
	mdl[@FBX_Model.vertex_buffer] = noone;
	mdl[@FBX_Model.material]    = -1;
	fbx_array_push(array,mdl);
}

function fbx_collect_deformer_data()
{
	/// @arg node
	/// @arg cluster_array
	/// @arg skin_array
	//!#import fbx_elem.* as FBX_Elem
	var node/*:FBX_Elem*/ = argument0, cluster_array = argument1, skin_array = argument2;
	var type = fbx_elem_prop_value(node, 2);
	switch type
	{
	    case "Cluster":
	        var dat/*:FBX_Deformer_Data*/ = array_create(FBX_Deformer_Data.sizeof);
	        dat[@FBX_Deformer_Data.uuid] = fbx_elem_uuid(node);
	        dat[@FBX_Deformer_Data.name] = fbx_elem_name(node);
	        dat[@FBX_Deformer_Data.type] = type;
	        dat[@FBX_Deformer_Data.indexes] = fbx_elem_property(node, "Indexes",0, [] );
	        dat[@FBX_Deformer_Data.weights] = fbx_elem_property(node, "Weights",0,[] );
	        dat[@FBX_Deformer_Data.transform] = fbx_array_from_list( fbx_elem_property(node, "Transform",0,[] ) );
	        dat[@FBX_Deformer_Data.transform_link] = fbx_array_from_list( fbx_elem_property(node, "TransformLink",0,[] ) );
	        // dat.transform_associate_model = fbx_array_from_list( node.property("TransformAssociateModel",0,[] ) );
	        fbx_array_push( cluster_array, dat );
	    break;
	    case "Skin":
	        var sdat/*:FBX_Skin_Data*/ = array_create(FBX_Skin_Data.sizeof);
	        sdat[@FBX_Skin_Data.uuid] = fbx_elem_uuid(node);
	        sdat[@FBX_Skin_Data.name] = fbx_elem_name(node);
	        sdat[@FBX_Skin_Data.type] = type;
	        sdat[@FBX_Skin_Data.deformers] = [];
	        fbx_array_push( skin_array, sdat);
	    break;
	}
}

function fbx_collect_anim_curve_node_data()
{
	/// @arg node
	/// @arg array
	//!#import fbx_elem.* as FBX_Elem
	var node/*:FBX_Elem*/ = argument0, array = argument1;
	var dat/*:FBX_AnimCurveNode_Data*/ = array_create(FBX_AnimCurveNode_Data.sizeof);
	dat[@FBX_AnimCurveNode_Data.uuid] = fbx_elem_uuid(node);
	dat[@FBX_AnimCurveNode_Data.name] = fbx_elem_name(node);
	dat[@FBX_AnimCurveNode_Data.type] = "";
	dat[@FBX_AnimCurveNode_Data.limb] = noone;
	dat[@FBX_AnimCurveNode_Data.transform] = [[],[],[]];
	fbx_array_push(array,dat);
}

function fbx_collect_anim_curve_data()
{
	/// @arg node
	/// @arg array
	//!#import fbx_elem.* as FBX_Elem
	var node/*:FBX_Elem*/ = argument0, array = argument1;
	var dat/*:FBX_AnimCurve_Data*/ = array_create(FBX_AnimCurve_Data.sizeof);
	dat[@FBX_AnimCurve_Data.uuid] = fbx_elem_uuid(node);
	dat[@FBX_AnimCurve_Data.name] = fbx_elem_name(node);
	var times_list = fbx_elem_property(node, "KeyTime",0,noone);
	if times_list == noone
	{
	    fbx_array_push(array,[[],[]]);
	    return false;
	}
	var size = ds_list_size(times_list);
	var ktimes = [];
	var i = -1;
	while ++i < size 
	{
	    var key = times_list[| i];
	    fbx_array_push(ktimes, key * FBX_KTIME);
	}
	dat[@FBX_AnimCurve_Data.keytime] = ktimes;
	dat[@FBX_AnimCurve_Data.keyvalue] = fbx_array_from_list( fbx_elem_property(node, "KeyValueFloat",0,[]) );
	fbx_array_push(array,dat);
}

function fbx_collect_alayer_data()
{
	/// @arg node
	/// @arg array
	//!#import fbx_elem.* as FBX_Elem
	var node/*:FBX_Elem*/ = argument0, array = argument1;
	var dat/*:FBX_ALayer_Data*/ = array_create(FBX_ALayer_Data.sizeof);
	dat[@FBX_ALayer_Data.uuid] = fbx_elem_uuid(node);
	dat[@FBX_ALayer_Data.name] = fbx_elem_name(node);
	dat[@FBX_ALayer_Data.nodes] = [];
	dat[@FBX_ALayer_Data.length] = 0;
	fbx_array_push(array,dat);
}

function fbx_collect_astack_data()
{
	/// @arg node
	/// @arg array
	//!#import fbx_elem.* as FBX_Elem
	var node/*:FBX_Elem*/ = argument0, array = argument1;
	var dat/*:FBX_AStack_Data*/ = array_create(FBX_AStack_Data.sizeof);
	dat[@FBX_AStack_Data.uuid] = fbx_elem_uuid(node);
	dat[@FBX_AStack_Data.name] = fbx_elem_name(node);
	fbx_array_push(array,dat);
}

function fbx_collect_node_attrib_data()
{
	/// @arg node
	/// @arg array
	//!#import fbx_elem.* as FBX_Elem
	var node/*:FBX_Elem*/ = argument0, array = argument1;
}

function fbx_collect_material_data()
{
	/// @arg node
	/// @arg array
	//!#import fbx_elem.* as FBX_Elem
	var node/*:FBX_Elem*/ = argument0, array = argument1;
}

function fbx_connections_map_create(fbx_conn)
{
	var conn_map = ds_map_create();
	var conn_map_reverse = ds_map_create();
	fbx_ds_list_add_map(FBX_TEMP_LISTS,conn_map);
	fbx_ds_list_add_map(FBX_TEMP_LISTS,conn_map_reverse);
	var elems = fbx_conn[| FBX_TREE];
	var i = -1;
	var size = ds_list_size(elems);
	while ++i < size
	{
	    var link = elems[| i];
	    var type = link[| FBX_TYPE];
	    var check_type = string_copy(type,2,2);
	    if check_type != "LL" then continue;
	    var data = link[| FBX_DATA];
	    var src = data[| 1];
	    var dst = data[| 2];
	    var src_list;
	    if ds_map_exists(conn_map,src) then src_list = conn_map[? src];
	    else {src_list = ds_list_create(); ds_map_add_list(conn_map,src,src_list);}
	    ds_list_add(src_list,[dst,link]);
	    var dst_list;
	    if ds_map_exists(conn_map_reverse,dst) then dst_list = conn_map_reverse[? dst];
	    else {dst_list = ds_list_create(); ds_map_add_list(conn_map_reverse,dst,dst_list);}
	    ds_list_add(dst_list,[src,link]);
	}
	return [conn_map,conn_map_reverse];
}

function fbx_collect_data()
{
	/// @arg fbx
	//!#import fbx_elem.* as FBX_Elem
	//!#import Blank
	var fbx/*:FBX_Elem*/ = argument0;
	;
	var nodes/*:FBX_Elem*/  = fbx_elem_find_first(fbx, "Objects");
	var fbx_connections/*:FBX_Elem*/ = fbx_elem_find_first(fbx, "Connections");
	if nodes == noone {show_error( "No Objects found", false);return noone;}
	var elems = fbx_elem_elems(nodes);
	if elems == noone {show_error( "No Elems found", false );return noone;}

	if fbx_connections == noone {show_error( "No Connections found", false);return noone;}
	var fbx_conns = fbx_connections_map_create(fbx_connections);
	var conn_map_rev = fbx_conns[1];

	var d/*:FBX_Pool*/;
	var i = -1;
	while ++i < FBX_Pool.sizeof
	{
	    d[i] = [];
	}

	d[@FBX_Pool.ties] = fbx_conns[1];
	d[@FBX_Pool.anim_index] = 0;

	var length = ds_list_size(elems);
	var i = -1;
	while ++i < length
	{
	    var node/*:FBX_Elem*/ = elems[|i];
	    var name = fbx_elem_id(node);
	    switch name 
	    {
	        case "Geometry":
	        fbx_collect_geometry_data(node,d[FBX_Pool.geom]);
	        break;
	        case "Model":
	        fbx_collect_model_data(node,d[FBX_Pool.limbs]);
	        break;
	        case "Pose":
	        fbx_collect_pose_data(node,d[FBX_Pool.pose]);
	        break;
	        case "Deformer":
	        fbx_collect_deformer_data(node,d[FBX_Pool.deformer],d[FBX_Pool.skin]);
	        break;
	        case "NodeAttribute":
	        fbx_collect_node_attrib_data(node,d[FBX_Pool.attrib]);
	        break;
	        case "Material":
	        fbx_collect_material_data(node,d[FBX_Pool.material]);
	        break;
	        case "AnimationStack":
	        fbx_collect_astack_data(node,d[FBX_Pool.astack]);
	        break;
	        case "AnimationLayer":
	        fbx_collect_alayer_data(node,d[FBX_Pool.alayer]);
	        break;
	        case "AnimationCurveNode":
	        fbx_collect_anim_curve_node_data(node,d[FBX_Pool.anode]);
	        break;
	        case "AnimationCurve":
	        fbx_collect_anim_curve_data( node, d[FBX_Pool.acurve] );
	        break;
	    }
	}
	return d;
}

function fbx_create_vertices_list(vertices, faces, indexes)
{
	var vert, size, iface, j, index, length, i, fid, px, py, pz;
	if faces == noone then return vertices;
	vert = fbx_ds_list_create();
	size = ds_list_size( faces );
	iface = [];
	j = -1;
	while ++j < size {
	    index = faces[|j];
	    if index < 0 {
	        index ^= -1;
	        fbx_array_push(iface,index);
	        length = array_length_1d(iface);
	        i = -1;
	        while ++i < length {
	            ds_list_add(indexes,iface[i]);
	            fid = iface[i] * 3;
	            px = vertices[|fid];
	            py = vertices[|fid+1];
	            pz = vertices[|fid+2];
	            ds_list_add(vert,px,py,pz);}
	        iface=[];}
	    else
	        fbx_array_push(iface,index);}

	return vert;
}
#endregion

#region fbx_connect.gml
function fbx_connect_curve_to_node()
{
	/// @arg data
	//!#import Blank
	var data/*:FBX_Connection_Data*/ = argument0;
	;
	var node/*:FBX_AnimCurveNode_Data*/ = data[FBX_Connection_Data.data];
	var link_elem = data[FBX_Connection_Data.link];
	var axis = fbx_elem_prop_value(link_elem,3);
	var src/*:FBX_AnimCurve_Data*/ = data[FBX_Connection_Data.res_data];
	var index = noone;
	switch axis
	{
	    case "d|X": index = 0;break;
	    case "d|Y": index = 1;break;
	    case "d|Z": index = 2;break;
	}
	if index >= 0
	{
	    var trn = node[FBX_AnimCurveNode_Data.transform];
	    trn[@index] = [src[FBX_AnimCurve_Data.keytime],src[FBX_AnimCurve_Data.keyvalue]];
	}
	return false; /// return false because we don't need to break while
}
function fbx_connect_limbs()
{
	/// @arg data
	//!#import FBX_Connection_Data
	var data/*:FBX_Connection_Data*/ = argument0;
	;
	var child/*:FBX_Model_Data*/ = data[FBX_Connection_Data.res_data];
	child[@FBX_Model_Data.parent] = data[FBX_Connection_Data.index];
	return false;
}
function fbx_connect_bones_to_deformers()
{
	/// @arg data
	//!#import Blank
	var data/*:FBX_Connection_Data*/ = argument0;
	;
	var deform/*:FBX_Deformer_Data*/ = data[FBX_Connection_Data.data];
	deform[@FBX_Deformer_Data.bone_index] = data[FBX_Connection_Data.res_index];
	var limb/*:FBX_Model_Data*/ = data[FBX_Connection_Data.res_data];
	limb[@FBX_Model_Data.deformer_index] = data[FBX_Connection_Data.index];
	return true; /// return false because we need to break while
}
function fbx_connect_clusters_to_skins()
{
	/// @arg data
	//!#import Blank
	var data/*:FBX_Connection_Data*/ = argument0;
	;
	var skin/*:FBX_Skin_Data*/ = data[FBX_Connection_Data.data];
	var child = data[FBX_Connection_Data.res_data];
	var sdefs = skin[FBX_Skin_Data.deformers];
	fbx_array_push(sdefs,child);
	return false; /// return false because we don't need to break while
}
function fbx_connect_anode_to_limbs()
{
/// @arg data
//!#import Blank
var data/*:FBX_Connection_Data*/ = argument0;
;
var link_elem = data[FBX_Connection_Data.link];
var type = fbx_elem_prop_value(link_elem,3);
var anode/*:FBX_AnimCurveNode_Data*/ = data[FBX_Connection_Data.res_data];
anode[@FBX_AnimCurveNode_Data.limb] = data[FBX_Connection_Data.index];
anode[@FBX_AnimCurveNode_Data.type] = type;
return false; /// return false because we don't need to break while
}

function fbx_connect_anode_to_alayer()
{
	/// @arg data
	//!#import FBX_Connection_Data
	var data/*:FBX_Connection_Data*/ = argument0;
	;
	var parent/*:FBX_ALayer_Data*/ = data[FBX_Connection_Data.data];
	var nodes = parent[FBX_ALayer_Data.nodes];
	var len = array_length_1d(nodes);
	var node/*:FBX_AnimCurveNode_Data*/ = data[FBX_Connection_Data.res_data];
	var index = node[FBX_AnimCurveNode_Data.limb];
	if index < 0 then return false;
	var ntrans;
	if index < len and index >= 0
	{
	    ntrans = nodes[index];
	    if not is_array(ntrans) 
	    {
	        ntrans = [[],[],[]];
	        nodes[@ index] = ntrans;
	    }
	}
	else
	{
	    ntrans = [[],[],[]];
	    nodes[@ index] = ntrans;
	}

	var j = noone;
	switch node[FBX_AnimCurveNode_Data.type]
	{
	    case "Lcl Translation":
	    j = 0;
	    break;
	    case "Lcl Rotation":
	    j = 1;
	    break;
	    case "Lcl Scaling":
	    j = 2;
	    break;
	}
	if j != noone
	{
	    ntrans[@ j] = node[FBX_AnimCurveNode_Data.transform];
	}

	//fbx_array_push(nodes,data.res_index);
	return false; /// return false because we don't need to break while
}

function fbx_connect_skin_to_geometry()
{
	/// @arg data
	//!#import FBX_Connection_Data
	var data/*:FBX_Connection_Data*/ = argument0;
	;
	var parent/*:FBX_Model*/ = data[FBX_Connection_Data.data];
	parent[@FBX_Model.skin] = data[FBX_Connection_Data.res_data];
	return false; /// return false because we don't need to break while
}

function fbx_connect_geom_to_limb()
{
	/// @arg data
	//!#import FBX_Connection_Data
	var data/*:FBX_Connection_Data*/ = argument0;
	;
	var geom/*:FBX_Model*/ = data[FBX_Connection_Data.res_data];
	geom[@FBX_Model.limb] = data[FBX_Connection_Data.index];
	return false; /// return false because we don't need to break while
}

function fbx_connect()
{
	/// @arg conn_type
	/// @arg a1
	/// @arg a2
	/// @arg conmap
	//!#import FBX_Connection_Data
	var conn_type = argument0, a1 = argument1, a2 = argument2, conmap = argument3;
	var size, i, chunk, uuid, con_list, list_size, j, con, child_uuid,
	    link_elem, conn_data/*:FBX_Connection_Data*/, con_res, res, script_index;
	script_index = asset_get_index(conn_type);
	size = array_length_1d(a1);
	i = -1;
	while ++i < size
	{
	    chunk = a1[i];
	    uuid = chunk[0];
	    con_list = conmap[? uuid];
	    if is_undefined(con_list) then continue;
	    list_size = ds_list_size(con_list);
	    j = -1;
	    while ++j < list_size
	    {
	        con = con_list[| j];
	        child_uuid  = con[0];
	        link_elem   = con[1];
	        res = fbx_data_by_uuid(a2,child_uuid);
	        if res == noone then continue;
	        conn_data = array_create(FBX_Connection_Data.sizeof);
	        conn_data[@FBX_Connection_Data.index]     = i;            // parent index
	        conn_data[@FBX_Connection_Data.data]      = chunk;        // parent itself
	        conn_data[@FBX_Connection_Data.link]      = link_elem;    // link elem
	        conn_data[@FBX_Connection_Data.res_index] = res[0];       // child index
	        conn_data[@FBX_Connection_Data.res_data]  = res[1];       // child itself
	        con_res = script_execute(script_index, conn_data);
	        if con_res then break;
	    }
	}
}

function fbx_data_by_uuid(array, uuid)
{
	var size = array_length_1d(array);
	var i = -1;
	while ++i < size
	{
	    var data = array[i];
	    if data[0] == uuid then return [i,data];
	}
	return noone;
}
#endregion

#region fbx_elem.gml
function fbx_elem_id(fbx_elem)
{
	///fbx_elem_id(fbx_elem)
	return argument0[| FBX_NAME];
}

function fbx_elem_type(fbx_elem)
{
	///fbx_elem_type(fbx_elem)
	return argument0[| FBX_TYPE];
}

function fbx_elem_props(fbx_elem)
{
	///fbx_elem_props(fbx_elem)
	return argument0[| FBX_DATA];
}

function fbx_elem_elems(fbx_elem)
{
	///fbx_elem_elems(fbx_elem)
	return argument0[| FBX_TREE];
}

function fbx_elem_uuid(fbx_elem)
{
	///fbx_elem_uuid(fbx_elem)
	var data = argument0[| FBX_DATA];
	return data[| 0];
}

function fbx_elem_name(fbx_elem)
{
///fbx_elem_name(fbx_elem)
	var data = argument0[| FBX_DATA];
	return fbx_string_to_array(data[| 1], chr(1) );
}

function fbx_elem_find_first(fbx_elem, key)
{
	var tree, length, i, elem, name;
	tree    = fbx_elem[|FBX_TREE];
	length  = ds_list_size(tree);
	i = -1;
	while ++i < length
	{
	    elem = tree[|i];
	    name = elem[|FBX_NAME];
	    if name == key then return elem;
	}
	return noone;
}

function fbx_elem_prop_first(elem, elem_id)
{
	if elem == noone then return noone;
	var elem_p70 = fbx_elem_find_first(elem,"Properties70");
	if elem_p70 == noone then return noone;
	var elems = elem_p70[| FBX_TREE];
	var i = -1;
	var size = ds_list_size( elems )
	while ++i < size
	{
	    var e = elems[| i];
	    var name = e[| FBX_NAME];
	    if name != "P" then continue;
	    var data = e[| FBX_DATA];
	    if data[| 0] == elem_id return data;
	}
	return noone;
}

function fbx_elem_prop_value(fbx_elem, i)
{
	var props = fbx_elem[| FBX_DATA];
	return props[| i];
}

function fbx_elem_prop_vector(fbx_elem, property, def_value)
{
	var res = fbx_elem_prop_first(fbx_elem,property);
	if res == noone then return def_value;
	return [res[| 4], res[| 5], res[| 6] ];
}

function fbx_elem_prop_bool(fbx_elem, property, def_value)
{
	var res = fbx_elem_prop_first(fbx_elem,property);
	if res == noone then return def_value;
	return res[| 4];
}

function fbx_elem_prop_enum(fbx_elem, property, def_value)
{
	var res = fbx_elem_prop_first(fbx_elem,property);
	if res == noone then return def_value;
	return res[| 4];
}

function fbx_elem_property(elem, child_name, prop_index, def_value)
{
	var e = fbx_elem_find_first(elem,child_name);
	if e == noone then return def_value
	return fbx_elem_prop_value(e,prop_index);
}
function fbx_elem_geom_layer()
{
	/// @arg node
	/// @arg layer_name
	/// @arg index_list
	//!#import fbx_elem.* as FBX_Elem
	var node/*:FBX_Elem*/ = argument0, layer_name = argument1, index_list = argument2;
	;
	var layr/*:FBX_Elem*/, map_inf, ref_inf, raw_name, ind_name, raw_list, ind_list,
	    type, list, direct_list;
	layr    = fbx_elem_find_first(node, layer_name );
	if layr == noone then return noone;
	map_inf = fbx_elem_property(layr, "MappingInformationType",0,"");
	ref_inf = fbx_elem_property(layr, "ReferenceInformationType",0,"");

	// trace( layer_name, map_inf, ref_inf );

	switch layer_name
	{
	    case "LayerElementNormal":
	        raw_name    = "Normals";
	        ind_name    = "NormalsIndex";
	        type        = 1;
	    break;
	    case "LayerElementBinormal":
	        raw_name    = "Binormals";
	        ind_name    = "BinormalsIndex";
	        type        = 1;
	    break;
	    case "LayerElementTangent":
	        raw_name    = "Tangents";
	        ind_name    = "TangentsIndex";
	        type        = 1;
	    break;
	    default:
	        raw_name    = "UV";
	        ind_name    = "UVIndex";
	        type        = 0;
	    break;
	}

	if ref_inf == "IndexToDirect"
	{
	    raw_list = fbx_elem_property(layr, raw_name,0,-4 );
	    ind_list = fbx_elem_property(layr, ind_name,0,-4 );
    
	    list = fbx_create_direct_list( raw_list, ind_list, type );
	}
	else // Direct
	{
	    list = fbx_elem_property(layr, raw_name,0,-4 );
	}

	if map_inf == "ByPolygonVertex"
	{
	    return list;
	}
	else // ByVertex or ByVertice
	{
	    direct_list = fbx_create_direct_list( list, index_list, type );
	    return direct_list;
	}
}

function fbx_create_direct_list(raw_list, ind_list, type)
{
	var list, length, i,j, uvid, px,py,pz;
	list = fbx_ds_list_create();
	length = ds_list_size(ind_list);
	i = -1;

	if type == 1
	{
	    while ++i < length 
	    {
	        uvid = ind_list[|i] * 3;
	        px = raw_list[|uvid];
	        py = raw_list[|uvid+1];
	        pz = raw_list[|uvid+2];
	        ds_list_add(list,px,py,pz);
	    }
	}
	else
	{
	    while ++i < length 
	    {
	        uvid = ind_list[|i] * 2;
	        uvid = max( uvid, 0);
	        px = raw_list[|uvid];
	        py = raw_list[|uvid+1];
	        ds_list_add(list,px,py);
	    }
	}

	return list;
}
#endregion

#region fbx_buffer.gml
function fbx_build_anim_buffer()
{
	var fbx/*:FBX_Pool*/ = argument0;
	if fbx[FBX_Pool.anim_count] == 0 return false;
	var layers, lay/*:FBX_ALayer_Data*/, anim, len, size, time, 
	    i, nodes, mtx, tick, buffer;
	layers = fbx[FBX_Pool.alayer];
	lay = layers[0];
	nodes = lay[FBX_ALayer_Data.nodes];
	size = array_length_1d(nodes);
	len = lay[FBX_ALayer_Data.length];
	tick = .015;
	buffer = buffer_create( 1, buffer_grow, 1 );
	buffer_write( buffer, buffer_u32, size );
	buffer_write( buffer, buffer_u32, floor( len / tick ) );
	i = -1;
	while ++i < size
	{
	    time = 0;
	    anim = nodes[i];
	    while time + tick < len
	    {
	        if not is_array(anim) then mtx = matrix_build_identity();
	        else mtx = fbx_frame_matrix( anim, time );
	        fbx_buffer_write_array( buffer, buffer_f32, mtx );
	        time += tick;
	    }
	}

	size = buffer_tell( buffer );
	return [buffer, size];
}

function fbx_vbuff_save()
{
	var fbx/*:FBX_Pool*/ = argument0;
	var name = argument1;
	var geoms = fbx[FBX_Pool.geom];

	var len = array_length_1d( geoms );
	var header_size = ( len * 2  + 1 ) * 4;
	var total_size = header_size;
	var geom/*:FBX_Model*/;
	var offsets, gbuff, size, offset, offset_pair;
	offsets = [];
	var i = -1;
	while ++i < len
	{
	    geom = geoms[i];
	    gbuff = geom[FBX_Model.buffer];
	    buffer_seek( gbuff, buffer_seek_end, 0 );
	    size = buffer_tell( gbuff );
	    fbx_array_push( offsets, [total_size, size] );
	    total_size += size;
	}

	var buffer = buffer_create( 1, buffer_grow, 1 );
	var buff_pair;
	buff_pair = fbx_build_bind_buffer(fbx);
	fbx_buffer_insert( buffer, FBX_BuffTypes.bind, buff_pair[1], buff_pair[0] );

	buff_pair = fbx_build_skin_buffer(fbx);
	fbx_buffer_insert( buffer, FBX_BuffTypes.skin, buff_pair[1], buff_pair[0] );

	buff_pair = fbx_build_anim_buffer(fbx);
	if buff_pair != false 
	{
	    fbx_buffer_insert( buffer, FBX_BuffTypes.anim, buff_pair[1], buff_pair[0] );
	}

	var i = -1;
	while ++i < len
	{
	    geom = geoms[i];
	    gbuff = geom[FBX_Model.buffer];
	    buffer_seek( gbuff, buffer_seek_end, 0 );
	    size = buffer_tell( gbuff );
	    fbx_buffer_insert( buffer, FBX_BuffTypes.geom, size, gbuff );
	}

	size = buffer_tell( buffer );
	var save_buff = buffer_create( size, buffer_fixed, 1 );
	buffer_copy( buffer, 0, size, save_buff, 0 );

	buffer_save( save_buff, name);
	buffer_delete(save_buff);
	buffer_delete(buffer);

	// var buff = buffer_create( total_size, buffer_fixed, 1);
	// buffer_write( buff, buffer_u32, len );
	// var i = -1;
	// while ++i < len
	// {
	//     geom = geoms[i];
	//     gbuff = geom.buffer;
	//     offset_pair = offsets[i];
	//     offset = offset_pair[0];
	//     size    = offset_pair[1];
	//     buffer_write( buff, buffer_u32, offset );
	//     buffer_write( buff, buffer_u32, size );
	//     buffer_copy( gbuff, 0, size, buff, offset );
	// }

	// buffer_save( buff, "testmeshes.msh");
	// buffer_delete(buff);
}
function fbx_buffer_write_array(buffer, type, array)
{
	var len, i;
	len = array_length_1d( array );
	i = -1;
	while ++i < len
	{
	    buffer_write( buffer, type, array[ i ] );
	}
	return true;
}
function fbx_buffer_insert(buffer, type, size, source)
{
	var offset;
	buffer_write( buffer, buffer_u32, size );
	buffer_write( buffer, buffer_u32, type );
	offset = buffer_tell( buffer );
	buffer_copy( source,0,size,buffer,offset);
	offset += size;
	buffer_seek( buffer, buffer_seek_start, offset );
	buffer_delete( source );
}
function fbx_build_skin_buffer()
{
	//!#import Blank
	var fbx/*:FBX_Pool*/ = argument0;
	var geoms, skin/*:FBX_Skin_Data*/, geom/*:FBX_Model*/, len, i, buffer, type,
	    deformers, deformer/*:FBX_Deformer_Data*/, j, size, index, bsize;
	geoms = fbx[FBX_Pool.geom];
	type = buffer_u8;
	buffer = buffer_create( 1, buffer_grow, 1 );
	len = array_length_1d( geoms );
	buffer_write( buffer, type, len );
	i = -1;
	while ++i < len
	{
	    geom = geoms[ i ];
	    skin = geom[FBX_Model.skin];
	    if skin == noone
	    {
	        buffer_write( buffer, type, 0 );
	        continue;
	    }
	    deformers = skin[FBX_Skin_Data.deformers];
	    size = array_length_1d( deformers );
	    buffer_write( buffer, type, size );
	    j = -1;
	    while ++j < size
	    {
	        deformer = deformers[ j ];
	        index = deformer[FBX_Deformer_Data.bone_index];
	        buffer_write( buffer, type, index );
	    }
	}
	bsize = buffer_tell( buffer );
	return [buffer, bsize];
}
function fbx_build_bind_buffer()
{
	//!#import Blank
	var fbx/*:FBX_Pool*/ = argument0;
	var limbs, limb/*:FBX_Model_Data*/, len, i, name, 
	    bind, parent, lbuff, lsize;
	lbuff = buffer_create( 1, buffer_grow, 1 );
	limbs = fbx[FBX_Pool.limbs];
	len = array_length_1d( limbs );
	buffer_write( lbuff, buffer_f32, len );
	i = -1;
	while ++i < len
	{
	    limb = limbs[i];
	    name = limb[FBX_Model_Data.name];
	    name = name[0];
	    bind = limb[FBX_Model_Data.bind_imt];
	    parent = limb[FBX_Model_Data.parent];
	    buffer_write( lbuff, buffer_f32, parent );
	    fbx_buffer_write_array( lbuff, buffer_f32, bind );
	}
	lsize = buffer_tell( lbuff );

	return [lbuff, lsize];
}
function fbx_object_anim_save()
{
	var fbx/*:FBX_Pool*/ = argument0, filename = argument1, obj_id = argument2;
	var limb/*:FBX_Model_Data*/;
	var limbs = fbx[FBX_Pool.limbs];
	var len = array_length_1d(limbs);

	var mtx = fbx[FBX_Pool.bone_mtx];
	var len = fbx[FBX_Pool.anim_length];
	var time = fbx[FBX_Pool.anim_time];
	var obj_mtx = [];
	var tick = 1/60;
	var frames = floor( len / tick );
	var framesize = 16 * 4;
	var bsize = framesize * frames + 4;
	var buff = buffer_create( bsize, buffer_fixed, 4 );
	buffer_write( buff, buffer_f32, frames );
	var px, py, pz, vx, vy, vz;
	var mat_rot = matrix_build(0,0,0,90,0,0,1,1,1);
	var pos, view;
	var i = -1;
	trace( "FRAMES", frames)
	while ++i < frames
	{
	    trace( i );
	    fbx_animate(fbx, tick);
	    obj_mtx = mtx[obj_id];
	    //obj_mtx = matrix_multiply( obj_mtx, mat_rot );
	    fbx_buffer_write_array( buff, buffer_f32, obj_mtx );
	    time = fbx[FBX_Pool.anim_time];
	}
	buffer_save( buff, filename );
	buffer_delete( buff );
	return true;
}
#endregion