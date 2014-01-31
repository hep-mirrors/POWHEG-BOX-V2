c -*- Fortran -*-

c Maximum number of weights
      integer lhrwgt_maxnheader
      parameter (lhrwgt_maxnheader=200)
      integer  lhrwgt_max_header_columns
      parameter ( lhrwgt_max_header_columns=200)
      character * (lhrwgt_max_header_columns)
     1     lhrwgt_header(lhrwgt_maxnheader)
      character * 100 lhrwgt_id,lhrwgt_group_name,
     1     lhrwgt_group_combine
      integer lhrwgt_nheader
      character * 400 lhrwgt_descr
      common/pwhg_lhrwgt/lhrwgt_nheader,lhrwgt_header,
     1     lhrwgt_id,lhrwgt_descr,lhrwgt_group_name,
     2     lhrwgt_group_combine
