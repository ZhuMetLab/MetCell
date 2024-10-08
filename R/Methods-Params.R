# project: MetCell
# File name: MethodsParams.R
# Created by: Yandong Yin
# Contact: yddream@gmail.com
# Created on: 2022/1/17 9:50
# Copyright (c) 2022- ZhuMSLab ALL right reserved

#' @export
setMethod("as.list", signature(x = "MetCellParam"), function(x, ...) {
  return(.param2list(x))
})

# The 'setAs' method.
setAs("MetCellParam" ,"list", function(from){
  return(.param2list(from))
})