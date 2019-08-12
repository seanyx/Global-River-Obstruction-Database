#' Reconstruct centerline (connectivity) from a table of coordinates of the centerline pixel
#'
#' @param cl_raw A tibble contains the UTM coordinates of the centerline.
#' @param reach_length A fixed reach length. Used to define simple reaches.
#' @return A tibble contains the segment id and reach id for the centerline.
#' @examples
#'

recon_cl = function(cl_raw, reach_length = 1000) {
  
  # convert to tibble
  cl_raw = as.tibble(cl_raw) %>%
    mutate(cl_id = 1:nrow(cl_raw))
  
  # apply kd-tree nearest neighbor search
  xyref = cl_raw %>% dplyr::select(x, y, cl_id)
  coords = xyref %>% dplyr::select(x, y)
  nearest = nn2(coords, coords)
  
  ind = as.tibble(nearest$nn.idx)
  dist = as.tibble(nearest$nn.dists)
  
  # construct nbTable (neighbor defined as distance smaller than 45 meters)
  # Var name|description
  # cl_id|unique centerline point index
  # nn_id|neighbor's cl_id
  # dist|distance in meters between points cl_id and nn_id
  
  t1 = Sys.time()
  
  print("Creating neighborhood table.")
  
  ## nn list cl_id and its number of neighbors within 45 meters
  nn = bind_cols(tibble(nn = apply(dist, 1, function(row) {return(sum(row < 45) - 1)})),
                 xyref %>% select(cl_id))
  
  nbTable = foreach(i = 1:nrow(nn), .combine = bind_rows) %do% {
    cl_id = i
    nni = nn$nn[i]
    if (nni != 0) {
      nn_id = (ind[i, 2:(nni + 1)] %>% gather())$value
      nn_dist = (dist[i, 2:(nni + 1)] %>% gather())$value
      tibble(cl_id = rep(cl_id, nni), nn_id = nn_id, dist = nn_dist)
    }
  }
  nbTableCost = Sys.time() - t1
  
  # define ends as only having exact one neighbor
  ends = nn %>% filter(nn == 1) %>% select(cl_id)
  # define joints as having more than 2 neighbors
  joints = nn %>% filter(nn > 2)
  
  # define new ends exposed after removing joints,
  # they are the neighbors of the joints
  newEnds = tibble(cl_id = setdiff((joints %>% left_join(nbTable, by = "cl_id"))$nn_id,
                                   joints$cl_id))
  
  # add the new ends
  ends = bind_rows(ends, newEnds)
  
  
  # remove rows of nbTable that have joints as the cl_id, so that each
  # 1. search beginning at one end points will not pass the joints,
  # 2. the resulting segments will always have two ends and no joints
  
  nbTable_jRemoved = nbTable %>% left_join(joints, by = "cl_id") %>%
    filter(is.na(nn)) %>%
    select(cl_id, nn_id, dist)
  
  # loop through ends, use nbTable_jRemoved to reconstruct connectivity
  t1 = Sys.time()
  segIdTable = tibble(cl_id = NA, segId = NA, segOrder = NA)
  segID = 1
  
  print("Reconstruct segments")
  
  for (i in 1:nrow(ends)) {
    
    s = ends[i, ]
    
    # check if this end point has been assigned segId,
    # if true, skip to next end point
    
    nbs = nbTable_jRemoved %>% right_join(s, by = "cl_id") %>% select(nn_id)
    if (is.na(nbs$nn_id)[1]) {
      print(paste0(i, 'th end point already assigned segID, so skipped.'))
      next
    }
    
    # if false, continue
    n1 = 1
    n2 = 2
    
    # while loop will add new neighbors of current segment incrementally under
    # the same segId, until there is no new neighbor to add
    while (n1 < n2) {
      n1 = nrow(s)
      nbs = nbTable_jRemoved %>% right_join(s, by = "cl_id") %>% select(nn_id)
      s = tibble(cl_id = c(s$cl_id, setdiff(nbs$nn_id, s$cl_id)))
      n2 = nrow(s)
    }
    ns = nrow(s)
    s = s[-ns, ] # remove the last s which is NA, s contains the cl_id for one segment
    
    # use unique combination of cl_id and nn_id to extract distance from nbTable, add
    # that to s
    s = bind_cols(s[1:(ns - 1), 1], s[2:ns, 1] %>% rename(nn_id = cl_id)) %>%
      left_join(nbTable, by = c("cl_id", "nn_id")) %>%
      select(cl_id, dist) %>%
      right_join(s, by = "cl_id")
    
    # add s to the segIdTable, increase segID by 1, add the order of the cl_id in the s
    segIdTable = bind_rows(segIdTable, s %>% mutate(segId = segID, segOrder = 1:nrow(s)))
    segID = segID + 1
    print(paste0(i, 'th point done with ', ns - 1, ' points.'))
    
    # remove the points that has been assigned reachId
    ids_toKeep = setdiff((nbTable_jRemoved$cl_id), s$cl_id)
    nbTable_jRemoved = nbTable_jRemoved %>% inner_join(tibble(cl_id = ids_toKeep), by = "cl_id")
  }
  t2 = Sys.time()
  print(t2 - t1) # 4min
  
  # remove the first row of the segIdTable
  segIdTable = segIdTable %>% filter(!is.na(cl_id))
  
  # add segment ID info to the clRaw input
  cl = xyref %>% right_join(segIdTable, by = "cl_id")
  
  # now cl has segment ID and oneEnd distance, lets assign reach ID to it
  # a quick method for doing so is to assign a unique reachID for every ~1km of the segment
  
  
  # assign reachId for each cl_id point in cl
  cl = cl %>%
    filter(!is.na(dist)) %>%
    group_by(segId) %>%
    mutate(ddist = cumsum(dist),
           reachId = segId * 1000 + ddist %/% reach_length) %>%
    ungroup() %>%
    left_join(cl_raw %>% select(lon, lat, cl_id), by = "cl_id")
  
  return(cl)
}


#' Reconstruct centerline (connectivity) from a table of coordinates of the centerline pixel
#'
#' @param cl_raw A tibble contains the UTM coordinates of the centerline.
#' @param reach_length A fixed reach length. Used to define simple reaches.
#' @return A tibble contains the segment id and reach id for the centerline.
#' @examples
#'

recon_cl_backup = function(cl_raw, reach_length = 1000) {
  
  # convert to tibble
  cl_raw = as.tibble(cl_raw) %>%
    mutate(cl_id = 1:nrow(cl_raw))
  
  # apply kd-tree nearest neighbor search
  xyref = cl_raw %>% dplyr::select(x, y, cl_id)
  coords = xyref %>% dplyr::select(x, y)
  nearest = nn2(coords, coords)
  
  ind = as.tibble(nearest$nn.idx)
  dist = as.tibble(nearest$nn.dists)
  
  # construct nbTable (neighbor defined as distance smaller than 45 meters)
  # Var name|description
  # cl_id|unique centerline point index
  # nn_id|neighbor's cl_id
  # dist|distance in meters between points cl_id and nn_id
  
  t1 = Sys.time()
  
  print("Creating neighborhood table.")
  
  ## nn list cl_id and its number of neighbors within 45 meters
  nn = bind_cols(tibble(nn = apply(dist, 1, function(row) {return(sum(row < 45) - 1)})),
                 xyref %>% select(cl_id))
  
  nbTable = foreach(i = 1:nrow(nn), .combine = bind_rows) %do% {
    cl_id = i
    nni = nn$nn[i]
    if (nni != 0) {
      nn_id = (ind[i, 2:(nni + 1)] %>% gather())$value
      nn_dist = (dist[i, 2:(nni + 1)] %>% gather())$value
      tibble(cl_id = rep(cl_id, nni), nn_id = nn_id, dist = nn_dist)
    }
  }
  nbTableCost = Sys.time() - t1
  
  # define ends as only having exact one neighbor
  ends = nn %>% filter(nn == 1) %>% select(cl_id)
  # define joints as having more than 2 neighbors
  joints = nn %>% filter(nn > 2)
  
  # define new ends exposed after removing joints,
  # they are the neighbors of the joints
  newEnds = tibble(cl_id = setdiff((joints %>% left_join(nbTable, by = "cl_id"))$nn_id,
                                   joints$cl_id))
  
  # add the new ends
  ends = bind_rows(ends, newEnds)
  
  
  # remove rows of nbTable that have joints as the cl_id, so that each
  # 1. search beginning at one end points will not pass the joints,
  # 2. the resulting segments will always have two ends and no joints
  
  nbTable_jRemoved = nbTable %>% left_join(joints, by = "cl_id") %>%
    filter(is.na(nn)) %>%
    select(cl_id, nn_id, dist)
  
  # loop through ends, use nbTable_jRemoved to reconstruct connectivity
  t1 = Sys.time()
  segIdTable = tibble(cl_id = NA, segId = NA, segOrder = NA)
  segID = 1
  
  print("Reconstruct segments")
  
  for (i in 1:nrow(ends)) {
    
    s = ends[i, ]
    
    # check if this end point has been assigned segId,
    # if true, skip to next end point
    
    nbs = nbTable_jRemoved %>% right_join(s, by = "cl_id") %>% select(nn_id)
    if (is.na(nbs$nn_id)[1]) {
      print(paste0(i, 'th end point already assigned segID, so skipped.'))
      next
    }
    
    # if false, continue
    n1 = 1
    n2 = 2
    
    # while loop will add new neighbors of current segment incrementally under
    # the same segId, until there is no new neighbor to add
    while (n1 < n2) {
      n1 = nrow(s)
      nbs = nbTable_jRemoved %>% right_join(s, by = "cl_id") %>% select(nn_id)
      s = tibble(cl_id = c(s$cl_id, setdiff(nbs$nn_id, s$cl_id)))
      n2 = nrow(s)
    }
    ns = nrow(s)
    s = s[-ns, ] # remove the last s which is NA, s contains the cl_id for one segment
    
    # use unique combination of cl_id and nn_id to extract distance from nbTable, add
    # that to s
    s = bind_cols(s[1:(ns - 1), 1], s[2:ns, 1] %>% rename(nn_id = cl_id)) %>%
      left_join(nbTable, by = c("cl_id", "nn_id")) %>%
      select(cl_id, dist) %>%
      right_join(s, by = "cl_id")
    
    # add s to the segIdTable, increase segID by 1, add the order of the cl_id in the s
    segIdTable = bind_rows(segIdTable, s %>% mutate(segId = segID, segOrder = 1:nrow(s)))
    segID = segID + 1
    print(paste0(i, 'th point done with ', ns - 1, ' points.'))
    
    # remove the points that has been assigned reachId
    ids_toKeep = setdiff((nbTable_jRemoved$cl_id), s$cl_id)
    nbTable_jRemoved = nbTable_jRemoved %>% inner_join(tibble(cl_id = ids_toKeep), by = "cl_id")
  }
  t2 = Sys.time()
  print(t2 - t1) # 4min
  
  # remove the first row of the segIdTable
  segIdTable = segIdTable %>% filter(!is.na(cl_id))
  
  # add segment ID info to the clRaw input
  cl = xyref %>% right_join(segIdTable, by = "cl_id")
  
  # now cl has segment ID and oneEnd distance, lets assign reach ID to it
  # a quick method for doing so is to assign a unique reachID for every ~1km of the segment
  
  
  # assign reachId for each cl_id point in cl
  cl = cl %>%
    filter(!is.na(dist)) %>%
    group_by(segId) %>%
    mutate(ddist = cumsum(dist),
           reachId = segId * 1000 + ddist %/% reach_length) %>%
    ungroup() %>%
    left_join(cl_raw %>% select(lon, lat, cl_id), by = "cl_id")
  
  return(cl)
}