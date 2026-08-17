%======================================================================
%                    E N D _ P R O C E S S I N G _ S T E P . M 
%                    doc: Fri Jun 25 16:17:17 2004
%                    dlm: Fri Jul 23 19:49:17 2004
%                    (c) 2004 ladcp@
%                    uE-Info: 25 46 NIL 0 0 72 2 2 8 NIL ofnI
%======================================================================

% finish processing step (in [process_cast.m])

if pcs.cur_step > 0
  disp(sprintf('==> STEP %d TOOK %.1f seconds',pcs.cur_step,toc-last_toc));
end
last_toc = toc;
if pcs.stop > 0 & pcs.cur_step >= pcs.target_begin_step
  if pcs.stop == 1, pcs.stop = 0; end
  disp(sprintf('entering DEBUG mode AFTER step %d (%s)',pcs.cur_step,pcs.step_name));
  disp(sprintf('(next stop = %d; type "return" to continue, "dbquit" to abort)',pcs.stop));
  keyboard;
  more off; % just in case...
end

if any(ismember(pcs.cur_step,p.checkpoints))
  disp(sprintf('SAVING CHECKPOINT %s_%d',f.checkpoints,pcs.cur_step));
  eval(sprintf('save %s_%d.mat',f.checkpoints,pcs.cur_step));
end

