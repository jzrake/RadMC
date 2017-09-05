


class Task(object):
    """A base class for tasks."""
    def __init__(self):
        pass

    def run(self, status, repetition):
        pass

    def get_recurrence(self):
        return Recurrence(0.0, 0.0, 1)

    def get_name(self):
        return self.__class__.__name__




class Recurrence(object):

    def __init__(self, phys_time_interval, wall_time_interval=0.0, iteration_interval=0):
        self.phys_time_interval = phys_time_interval
        self.wall_time_interval = wall_time_interval
        self.iteration_interval = iteration_interval

        self.next_phys_time = 0.0
        self.next_wall_time = 0.0
        self.next_iteration = 0
        self.repetition = 0
        self.skip_next = False

    def is_due(self, status):
        if self.phys_time_interval > 0.0 and self.next_phys_time <= status.simulation_time: return 'P'
        if self.wall_time_interval > 0.0 and self.next_wall_time <= status.wall_minutes: return 'W'
        if self.iteration_interval > 0   and self.next_iteration <= status.simulation_iter: return 'I'
        return None

    def update(self, reason):
        self.repetition += 1
        self.skip_next = False
        if reason == 'P': self.next_phys_time += self.phys_time_interval
        if reason == 'W': self.next_wall_time += self.wall_time_interval
        if reason == 'I': self.next_iteration += self.iteration_interval




class TaskDescription(object):

    def __init__(self, task, recurrence, name):
        self.task = task
        self.recurrence = recurrence
        self.name = name




class TaskScheduler(object):
    """A class for dispatching tasks within a main loop."""
    def __init__(self):
        self.tasks = dict()

    def schedule(self, task):
        name = task.get_name()
        if name in self.tasks:
            raise KeyError(name + " already scheduled")
        self.tasks[name] = TaskDescription(task, task.get_recurrence(), name)

    def skip_next(self, name):
        self.tasks[name].recurrence.skip_next = True
    
    def dispatch(self, status):
        for name, task_description in self.tasks.iteritems():

            recr = task_description.recurrence
            task = task_description.task

            reason = recr.is_due(status)

            if reason:
                if not recr.skip_next:
                    task.run(status, recr.repetition)
                recr.update(reason)




class Status(object):
    """Status class holding iteration, simulation time, and wall time."""
    import time
    def __init__(self):
        self.simulation_iter = 0
        self.simulation_time = 0.0
        self.wall_minutes = 0.0
        self.time_started = self.time.time()

    def step(self, dt):
        self.simulation_time += dt
        self.simulation_iter += 1
        self.wall_minutes = (self.time.time() - self.time_started) / 60




if __name__ == "__main__":
    class MyTask(Task):
        def run(self, status, repetition):
            print 'running task', repetition, status.simulation_iter
        def get_recurrence(self):
            return Recurrence(0.0, 0.0, 2)


    status = Status()
    scheduler = TaskScheduler()
    scheduler.schedule(MyTask())

    for i in range(10):
        scheduler.dispatch(status)
        status.step(0.1)
